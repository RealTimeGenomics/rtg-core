/*
 * Copyright (c) 2014. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
 */
package com.rtg.variant.eval;

import java.util.EnumSet;

import com.rtg.mode.DnaUtils;
import com.rtg.util.Utils;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.intervals.SequenceNameLocus;
import com.rtg.util.intervals.SequenceNameLocusComparator;
import com.rtg.vcf.VcfRecord;
import com.rtg.vcf.VcfUtils;


/**
 * Class for holding information parsed from a SNP detection file
 */
public class DetectedVariant implements Comparable<DetectedVariant>, Variant {

  static final SequenceNameLocusComparator NATURAL_COMPARATOR = new SequenceNameLocusComparator();

  private final String mSequenceName;
  private final double mSortValue;
  private final int mStart;
  private final int mEnd;
  private final byte[][] mPrediction;
  private final boolean mPhased;
  private final EnumSet<RocFilter> mFilters = EnumSet.noneOf(RocFilter.class);

  /**
   * Construct a Detected Variant by inspecting a <code>VcfRecord</code> object. A string is supplied as the output representation of this object
   * @param rec VCF record to process
   * @param sampleNo the sample column number (starting from 0) for multiple sample variant calls
   * @param extractor ROC value extractor implementation to use
   * @param squashPloidy if true, convert heterozygous calls to homozygous alternative.
   */
  public DetectedVariant(VcfRecord rec, int sampleNo, RocSortValueExtractor extractor, boolean squashPloidy) {
    mSequenceName = new String(rec.getSequenceName().toCharArray()); // Prevent entire VCF record string from being kept
    mSortValue = extractor.getSortValue(rec, sampleNo);
    final boolean hasPreviousNt = VcfUtils.hasRedundantFirstNucleotide(rec);
    mStart = rec.getStart() + (hasPreviousNt ? 1 : 0);
    mEnd = rec.getEnd();
    final String gt = rec.getFormatAndSample().get(VcfUtils.FORMAT_GENOTYPE).get(sampleNo);
    final int[] gtArray = VcfUtils.splitGt(gt);
    mPhased = gt.contains("" + VcfUtils.PHASED_SEPARATOR);
    assert gtArray.length == 1 || gtArray.length == 2; //can only handle haploid or diploid
    mPrediction = new byte[VcfUtils.isHomozygous(rec, sampleNo) || squashPloidy ? 1 : 2][];
    if (squashPloidy) {
      final int gtId = gtArray.length == 1 ? gtArray[0] : Math.max(gtArray[0], gtArray[1]);
      mPrediction[0] = DnaUtils.encodeString(getAllele(rec, gtId, hasPreviousNt));
    } else {
      for (int i = 0; i < mPrediction.length; i++) {
        mPrediction[i] = DnaUtils.encodeString(getAllele(rec, gtArray[i], hasPreviousNt));
        //System.err.println(rec.getPosition() + " " + gtArray[i] + " " + getAllele(gtArray[i]));
      }
    }
    for (final RocFilter filter : RocFilter.values()) {
      if (filter.accept(rec, sampleNo)) {
        mFilters.add(filter);
      }
    }
  }

  private static String getAllele(VcfRecord rec, int allele, boolean hasPreviousNt) {
    if (allele == -1) {
      return "N"; // Missing allele
    }
    if (allele == 0) {
      if (hasPreviousNt) {
        if (rec.getRefCall().length() == 1) {
          return "";
        }
        return rec.getRefCall().substring(1);
      } else {
        return rec.getRefCall();
      }
    } else {
      if (rec.getAltCalls().size() < allele) {
        throw new NoTalkbackSlimException("Invalid allele number in record " + rec.toString());
      } else {
        final String localAllele = rec.getAltCalls().get(allele - 1);
        if (hasPreviousNt) {
          if (localAllele.length() == 1) {
            return "";
          }
          return localAllele.substring(1);
        } else {
          return localAllele;
        }
      }
    }
  }

  @Override
  public int compareTo(final DetectedVariant o2) {
    return NATURAL_COMPARATOR.compare(this, o2);
  }

  @Override
  public boolean equals(final Object o2) {
    if (this == o2) {
      return true;
    }
    if (o2 == null) {
      return false;
    }
    if (!(o2 instanceof DetectedVariant)) {
      return false;
    }
    return NATURAL_COMPARATOR.compare(this, (DetectedVariant) o2) == 0;
  }

  @Override
  public int hashCode() {
    return Utils.pairHash(getStart(), getSequenceName().hashCode());
  }

  @Override
  public String getSequenceName() {
    return mSequenceName;
  }

  @Override
  public String toString() {
    final StringBuilder sb = new StringBuilder();
    sb.append(getSequenceName()).append(":").append(getStart() + 1).append("-").append(getEnd() + 1).append(" ");
    sb.append(DnaUtils.bytesToSequenceIncCG(ntAlleleA()));
    final byte[] bAllele = ntAlleleB();
    if (bAllele != null) {
      sb.append(":").append(DnaUtils.bytesToSequenceIncCG(bAllele));
    }
    return sb.toString();
  }

  @Override
  public int getStart() {
    return mStart;
  }

  @Override
  public int getEnd() {
    return mEnd;
  }

  @Override
  public boolean overlaps(SequenceNameLocus other) {
    throw new UnsupportedOperationException();
  }

  @Override
  public boolean contains(String sequence, int pos) {
    throw new UnsupportedOperationException();
  }

  @Override
  public int getLength() {
    return getEnd() - getStart();
  }

  @Override
  public byte[] nt(boolean alleleA) {
    if (mPrediction.length == 2) {
      return mPrediction[alleleA ? 0 : 1];
    } else if (alleleA) {
      return mPrediction[0];
    } else {
      return null;
    }
  }

  @Override
  public byte[] ntAlleleA() {
    return nt(true);
  }

  @Override
  public byte[] ntAlleleB() {
    return nt(false);
  }

  @Override
  public boolean isPhased() {
    return mPhased;
  }

  /**
   * Get if the detected variant is accepted by the given <code>RocFilter</code>
   * @param filter the filter to check
   * @return true if the variant is accepted by the given filter, false otherwise
   */
  public boolean filterAccept(RocFilter filter) {
    return mFilters.contains(filter);
  }

  /**
   * Get the sort value for ROC curves
   * @return the sort value
   */
  public double getSortValue() {
    return mSortValue;
  }
}
