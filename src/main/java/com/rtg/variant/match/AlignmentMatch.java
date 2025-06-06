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
package com.rtg.variant.match;

import com.rtg.reader.Arm;
import com.rtg.reader.FastaUtils;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.Integrity;
import com.rtg.variant.MachineErrorChooserInterface;
import com.rtg.variant.PhredScaler;
import com.rtg.variant.VariantAlignmentRecord;
import com.rtg.mode.DNARange;
import com.rtg.mode.DNARangeNAT;
import com.rtg.variant.util.VariantUtils;

/**
 * Represent a match between a read (alignment) and reference.
 */
public class AlignmentMatch extends Match implements Integrity {

  private static final DNARange DNA_RANGE = DNARange.RANGE;

  private final double mMapError;
  private final int[] mReadNt;
  private final double[] mBaseError;
  private final double mQualityDefault;
  private final String mReadString;
  private final boolean mFixedLeft;
  private final boolean mFixedRight;
  private int mBasesLeftOfMatch;
  private int mBasesRightOfMatch;
  private int mSoftClipLeft;
  private int mSoftClipRight;
  private final VariantAlignmentRecord mAlignmentRecord;


  /**
   * For testing only, as it uses ascii qualities and does not use a <code>MachineErrorChooser</code>, so does no quality curve correction.
   * A new single nucleotide difference match.
   *
   * @param var associated alignment record (if any).
   * @param readNt read nucleotides
   * @param qScore Phred quality of read nucleotides
   * @param qDefault default Phred quality of read nucleotides Used if <code>qScore</code> null)
   * @param start zero-based offset into read nucleotides
   * @param length length of read region
   * @param readScore Phred quality of read mapping
   */
  public AlignmentMatch(final VariantAlignmentRecord var, final String readNt, final String qScore, final int qDefault, final int start, final int length, final int readScore) {
    this(var, null, readNt, FastaUtils.asciiToRawQuality(qScore), qDefault, start, length, readScore, true, true);
  }

  /**
   * A new single nucleotide difference match
   * @param var associated alignment record (if any).
   * @param chooser Machine error chooser (if null, no quality correction is done)
   * @param readNt read nucleotides
   * @param qScore binary Phred quality of read nucleotides
   * @param qDefault default Phred quality of read nucleotides (used if <code>qScore</code> null)
   * @param start zero-based offset into read nucleotides
   * @param length length of read region
   * @param readScore Phred quality of read mapping
   * @param fixedLeft true if read is fixed (anchored) at left end.
   * @param fixedRight true if read is fixed (anchored) at right end.
   */
  public AlignmentMatch(final VariantAlignmentRecord var, final MachineErrorChooserInterface chooser, final String readNt, final byte[] qScore, final int qDefault, final int start, final int length, final int readScore, final boolean fixedLeft, final boolean fixedRight) {
    super();
    mAlignmentRecord = var;
    final PhredScaler me = chooser == null ? null : chooser.machineErrors(var.getReadGroup(), var.isReadPaired());
    mMapError = VariantUtils.phredToProb(readScore);
    if (qScore == null) {
      mBaseError = null;
      mQualityDefault = VariantUtils.phredToProb(me == null ? qDefault : me.getScaledPhred((byte) qDefault, 0, Arm.LEFT));
    } else {
      mBaseError = new double[length];
      int l = length;
      double c = 0;
      if (l == 0) { // Compute average of the bounding bases
        if (start > 0) {
          c += VariantUtils.phredToProb(qScore[start - 1]);
          ++l;
        }
        if (start + length + 1 < qScore.length) {
          c += VariantUtils.phredToProb(qScore[start + length]);
          ++l;
        }
        if (l == 0) { // Fall back to default
          c += VariantUtils.phredToProb(me == null ? qDefault : me.getScaledPhred((byte) qDefault, 0, Arm.LEFT));
          ++l;
        }
      }
      for (int k = 0; k < mBaseError.length; ++k) {
        final int phred = qScore[k + start];
        assert 0 <= phred;
        mBaseError[k] = VariantUtils.phredToProb(phred);
        c += mBaseError[k];
      }
      mQualityDefault = c / l;
    }
    mFixedLeft = fixedLeft;
    mFixedRight = fixedRight;
    mReadNt = new int[length];
    final StringBuilder sb = new StringBuilder();
    for (int k = 0; k < length; ++k) {
      final int v = DNA_RANGE.valueOf(readNt.charAt(k + start));
      mReadNt[k] = v;
      sb.append(DNA_RANGE.toString(v));
    }
    mReadString = sb.toString();
  }

  /**
   * Get the associated SAM record.
   * @return the associated SAM record (may be null for legacy reasons).
   */
  public VariantAlignmentRecord alignmentRecord() {
    return mAlignmentRecord;
  }

  @Override
  public double baseError() {
    return mQualityDefault;
  }

  @Override
  public double baseError(final int index) {
    if (mBaseError == null) {
      return mQualityDefault;
    }
    return mBaseError[index];
  }

  @Override
  public int read(final int index) {
    final int res = mReadNt[index] - 1;
    assert DNARangeNAT.DNA.valid(res);
    return res;
  }

  public int getBasesLeftOfMatch() {
    return mBasesLeftOfMatch;
  }

  public void setBasesLeftOfMatch(int basesLeftOfMatch) {
    mBasesLeftOfMatch = basesLeftOfMatch;
  }

  public int getBasesRightOfMatch() {
    return mBasesRightOfMatch;
  }

  public void setBasesRightOfMatch(int basesRightOfMatch) {
    mBasesRightOfMatch = basesRightOfMatch;
  }

  /**
   * @return Position of first non soft clip base on read (0 based inclusive)
   */
  public int getSoftClipLeft() {
    return mSoftClipLeft;
  }

  /**
   * @param softClipLeft Position of first non soft clip base on read (0 based inclusive)
   */
  public void setSoftClipLeft(int softClipLeft) {
    mSoftClipLeft = softClipLeft;
  }

  /**
   * @return Position after last non soft clip base on read (i.e. 0 based exclusive)
   */
  public int getSoftClipRight() {
    return mSoftClipRight;
  }

  /**
   * @param softClipRight Position after last non soft clip base on read (i.e. 0 based exclusive)
   */
  public void setSoftClipRight(int softClipRight) {
    mSoftClipRight = softClipRight;
  }

  @Override
  public double mapError() {
    return mMapError;
  }

  @Override
  public int length() {
    return mReadNt.length;
  }

  @Override
  public String readString() {
    return mReadString;
  }

  @Override
  public boolean isFixedLeft() {
    return mFixedLeft;
  }

  @Override
  public boolean isFixedRight() {
    return mFixedRight;
  }

  @Override
  public boolean integrity() {
    Exam.assertEquals(mReadString, mReadString.length(), mReadNt.length);
    if (mBaseError != null) {
      Exam.assertEquals(mReadNt.length, mBaseError.length);
      Exam.assertProbabilities(mBaseError);
    }
    Exam.assertProbability(mMapError);
    Exam.assertProbability(mQualityDefault);
    return true;
  }

  @Override
  public boolean globalIntegrity() {
    return true;
  }
}
