/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
    for (int i = 0; helper.result() == 0 && i < getCoverageBitSet().length(); ++i) {
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
      for (int i = 0; i < cigar.length(); ++i) {
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
