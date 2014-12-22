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
package com.rtg.visualization;

import com.rtg.mode.DnaUtils;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.vcf.VcfRecord;
import com.rtg.vcf.VcfUtils;


/**
 * Class for holding information parsed from a SNP detection file
 */
public class AviewVariant {

  private final int mReferenceLength;
  private final int mPosition;
  private final boolean mIsFiltered;
  private final byte[][] mPrediction;

  /**
   * Construct a Detected Variant by inspecting a <code>VcfRecord</code> object. A string is supplied as the output representation of this object
   * @param rec VCF record to process
   * @param sampleNo the sample column number (starting from 0) for multiple sample variant calls
   */
  public AviewVariant(VcfRecord rec, int sampleNo) {
    final String gt = rec.getFormatAndSample().get(VcfUtils.FORMAT_GENOTYPE).get(sampleNo);
    final int[] gtArray = VcfUtils.splitGt(gt);
    if (gtArray.length > 2) {
      throw new NoTalkbackSlimException("Cannot handle calls with ploidy more than diploid.");
    }

    final boolean hasPreviousNt = VcfUtils.hasRedundantFirstNucleotide(rec);
    final String[] alleleStrings = VcfUtils.getAlleleStrings(rec, hasPreviousNt);
    mPrediction = new byte[VcfUtils.isHomozygous(rec, sampleNo) ? 1 : 2][];
    for (int i = 0; i < mPrediction.length; i++) {
      final int alleleId = gtArray[i];
      mPrediction[i] = DnaUtils.encodeString(alleleId == -1 ? "N" : alleleStrings[alleleId]);
    }
    mIsFiltered = rec.isFiltered();
    mReferenceLength = rec.getRefCall().length() - (hasPreviousNt ? 1 : 0);
    mPosition = rec.getOneBasedStart() + (hasPreviousNt ? 1 : 0);
  }


  /**
   * @return the position (1 based).
   */
  public int getPosition() {
    return mPosition;
  }

  /**
   * @return the length of the reference section
   */
  public int referenceLength() {
    return mReferenceLength;
  }

  /**
   * One allele of the variant as determined by allele parameter.
   * @param alleleA if true select the A allele.
   * @return the allele (may be null or zero length)
   */
  public byte[] nt(boolean alleleA) {
    if (mPrediction.length == 2) {
      return mPrediction[alleleA ? 0 : 1];
    } else if (alleleA) {
      return mPrediction[0];
    } else {
      return null;
    }
  }

  /**
   * The first allele of the variant. Always non null but may be zero length.
   * @return the A allele
   */
  public byte[] ntAlleleA() {
    return nt(true);
  }

  /**
   * Second allele of the variant. Will be null if homozygous and may be zero length.
   * @return the B allele, or null if homozygous
   */
  public byte[] ntAlleleB() {
    return nt(false);
  }

  /**
   * @return true if original record was filtered, false otherwise
   */
  public boolean isFiltered() {
    return mIsFiltered;
  }
}
