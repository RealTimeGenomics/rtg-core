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

package com.rtg.variant.sv;

import com.rtg.util.Utils;
import com.rtg.util.integrity.Exam;

/**
 * Square distribution, <code>rate1</code> at left and <code>rate2</code> at right.
 * <img src="doc-files/notch.jpg" alt="image">
 */
public final class DistributionStep extends Distribution {

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
