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
package com.rtg.variant.coverage;

import java.util.BitSet;

import com.rtg.sam.SamUtils;
import com.rtg.util.CompareHelper;

import htsjdk.samtools.SAMRecord;

/**
 * Hold the SAM record information we care about for the purposes of coverage calculations.
 */
public class CoverageReaderRecord extends AbstractMateInfoReaderRecord<CoverageReaderRecord> {

  private final BitSet mCoverageBitSet;
  private final int mIH;
  private final double mCoverageMultiplier;

  /**
   * @param sam record to convert
   * @param genome genome record applies to
   * @param includeDeletions if true, count deletions
   */
  public CoverageReaderRecord(SAMRecord sam, int genome, boolean includeDeletions) {
    super(sam, genome);
    final Integer ih = sam.getIntegerAttribute(SamUtils.ATTRIBUTE_IH);
    mIH = ih == null ? 1 : ih;
    mCoverageMultiplier = ih == null ? 1.0 : 1.0 / ih;
    mCoverageBitSet = parseCigar(sam, includeDeletions);
  }

  @Override
  public int compareTo(CoverageReaderRecord o) {
    return new CompareHelper().compare(getSequenceId(), o.getSequenceId())
          .compare(getStart(), o.getStart())
          .compare(getEnd(), o.getEnd())
          .compare(mGenome, o.mGenome)
          .compare(mMateSequenceId, o.mMateSequenceId)
          .compare(mFragmentLength, o.mFragmentLength)
          .compare(System.identityHashCode(this), System.identityHashCode(o)).result();
  }

  @Override
  public boolean equals(Object o) {
    if (o == null) {
      return false;
    }
    if (!(o instanceof CoverageReaderRecord)) {
      return false;
    }
    return compareTo((CoverageReaderRecord) o) == 0;
  }

  @Override
  public int hashCode() {
    return super.hashCode();
  }

  public int getIH() {
    return mIH;
  }

  public double getCoverageMultiplier() {
    return mCoverageMultiplier;
  }

  /**
   * @return a <code>BitSet</code> describing the coverage over the length of this record.
   */
  public BitSet getCoverageBitSet() {
    return mCoverageBitSet;
  }

  @Override
  public int disambiguateDuplicate(CoverageReaderRecord rec) {
    final CompareHelper helper = new CompareHelper().compare(getStart(), rec.getStart()).compare(getCoverageBitSet().length(), rec.getCoverageBitSet().length());
    for (int i = 0; helper.result() == 0 && i < getCoverageBitSet().length(); i++) {
      helper.compare(getCoverageBitSet().get(i), rec.getCoverageBitSet().get(i));
    }
    return helper.result();
  }

  private BitSet parseCigar(SAMRecord rec, boolean includeDeletions) {
    final BitSet bs = new BitSet();
    final String cigar = rec.getCigarString();
    if (!SamUtils.NO_CIGAR.equals(cigar)) {
      int tPos = 0;
      int n = 0;
      for (int i = 0; i < cigar.length(); i++) {
        final char c = cigar.charAt(i);
        if (Character.isDigit(c)) {
          n = 10 * n + c - '0';
        } else {
          switch (c) {
            case SamUtils.CIGAR_SAME_OR_MISMATCH:
            case SamUtils.CIGAR_SAME:
            case SamUtils.CIGAR_MISMATCH:
              // it is a match or mismatch, so increment our coverage counts
              bs.set(tPos, tPos + n, true);
              tPos += n;
              break;

            case SamUtils.CIGAR_GAP_IN_READ:
              // skip this region in the reference genome
              tPos += n;
              break;

            case SamUtils.CIGAR_DELETION_FROM_REF: // we record the delete at the read position just after the delete
              if (includeDeletions) {
                bs.set(tPos, tPos + n, true);
              }
              tPos += n;
              break;

            case SamUtils.CIGAR_INSERTION_INTO_REF:
            case SamUtils.CIGAR_SOFT_CLIP:
              // skip over these soft-clipped read bases.
            case SamUtils.CIGAR_HARD_CLIP:
              // unlike soft clipping, hard clipping does not appear in the read, so we just skip it.
            case SamUtils.CIGAR_PADDING:
              // padding is just to get multiple inserts aligned, so it does not increment either position.
              break;

            default:
              throw new RuntimeException("Unknown cigar code=" + c + " in record: " + rec.getSAMString());
          }
          n = 0;
        }
      }
    }
    return bs;
  }
}
