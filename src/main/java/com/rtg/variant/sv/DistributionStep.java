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

package com.rtg.variant.sv;

import com.rtg.util.Utils;
import com.rtg.util.integrity.Exam;

/**
 * Square distribution, <code>rate1</code> at left and <code>rate2</code> at right.
 * <img src="doc-files/notch.jpg" alt="image">
 */
public class DistributionStep extends Distribution {

  private final double mRate1;

  private final double mRate2;

  private final int mOffset;

  /**
   * @param lo low index of the distribution (inclusive).
   * @param hi high index of the distribution (exclusive).
   * @param offset in the window (lo to hi-1 inclusive)
   * @param rate1 rate for first (left) region (inclusive of offset).
   * @param rate2 rate for second (right) region.
   * @param reverse if true, swap <code>rate1</code> and <code>rate2</code>
   */
  public DistributionStep(final int lo, final int hi, final int offset, final double rate1, final double rate2, boolean reverse) {
    this(lo, hi, offset, reverse ? rate2 : rate1, reverse ? rate1 : rate2);
  }

  /**
   * @param lo low index of the distribution (inclusive).
   * @param hi high index of the distribution (exclusive).
   * @param offset in the window (-radius to radius -1 inclusive)
   * @param rate1 rate for first (left) region (inclusive of offset).
   * @param rate2 rate for second (right) region.
   */
  public DistributionStep(final int lo, final int hi, final int offset, final double rate1, final double rate2) {
    super(lo, hi);
    mOffset = offset;
    mRate1 = rate1;
    mRate2 = rate2;
    assert globalIntegrity();
  }

  @Override
  public Signal getSignalLn(SamCounts counts, String label) {
    return new SignalConstantLn(counts, this, label);
  }

  @Override
  protected double getValue(int index) {
    if (index <= mOffset) {
      return mRate1;
    }
    return mRate2;
  }

  /**
   * Returns the step offset
   * @return the step offset
   */
  public int getOffset() {
    return mOffset;
  }
  /**
   * Returns the first region rate
   * @return the first region rate
   */
  public double getRate1() {
    return mRate1;
  }
  /**
   * Returns the second region rate
   * @return the second region rate
   */
  public double getRate2() {
    return mRate2;
  }

  @Override
  public void toString(StringBuilder sb) {
    sb.append("Step:radius=").append(mHi).append(" offset=").append(mOffset).append(" rate1=").append(Utils.realFormat(mRate1, 4)).append(" rate2=").append(Utils.realFormat(mRate2, 4));
  }

  @Override
  public boolean integrity() {
    super.integrity();
    Exam.assertTrue(" " + mRate1, mRate1 >= 0.0 && Double.isFinite(mRate1));
    Exam.assertTrue(" " + mRate2, mRate2 >= 0.0 && Double.isFinite(mRate2));
    Exam.assertTrue(" offset=" + mOffset + " length=" + mHi, -mHi <= mOffset && mOffset < mHi);
    return true;
  }

}
