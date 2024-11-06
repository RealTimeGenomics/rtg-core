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

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.util.StringUtils;
import com.rtg.util.Utils;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 */
@TestClass(value = {"com.rtg.variant.sv.DistributionArrayTest"})
public abstract class Distribution extends IntegralAbstract {

  protected final int mLo;

  protected final int mHi;

  /**
   * @param lo low index of the distribution (inclusive).
   * @param hi high index of the distribution (exclusive).
   */
  public Distribution(int lo, int hi) {
    super();
    mLo = lo;
    mHi = hi;
    //System.err.println("lo=" + lo + " hi=" + hi);
  }

  /**
   * Get the value in the distribution.
   * @param index position in the distribution (from -radius to radius-1 inclusive).
   * @return the value.
   */
  public double get(final int index) {
    if (index < mLo || index >= mHi) {
      throw new IndexOutOfBoundsException("index=" + index + " lo=" + mLo + " hi=" + mHi);
    }
    return getValue(index);
  }

  /**
   * Return a log-based signal for this distribution with the given sam counts
   * @param counts the sam counts
   * @param label a textual label
   * @return the signal
   */
  public Signal getSignalLn(SamCounts counts, String label) {
    return new SignalDistributionLn(counts, this, label);
  }

  /**
   * Get a value in the distribution.
   * This version is local and does not check the index values.
   * @param index position in the distribution (0 based).
   * @return the value.
   */
  protected abstract double getValue(int index);

  /**
   * Lowest index in the distribution (inclusive).
   * @return the low index.
   */
  int lo() {
    return mLo;
  }

  /**
   * Highest index in the distribution (exclusive).
   * @return the high index.
   */
  int hi() {
    return mHi;
  }

  String dump() {
    final StringBuilder sb = new StringBuilder();
    for (int i = lo(); i < hi(); ++i) {
      sb.append(i).append(" ").append(Utils.realFormat(get(i), 4)).append(StringUtils.LS);
    }
    return sb.toString();
  }

  @Override
  public boolean integrity() {
    Exam.assertTrue(mLo < 0);
    Exam.assertTrue(mHi > 0);
    return true;
  }

  @Override
  public boolean globalIntegrity() {
    integrity();
    for (int i = lo(); i < hi(); ++i) {
      final double v = getValue(i);
      Exam.assertTrue(Double.isFinite(v));
    }
    return true;
  }

}
