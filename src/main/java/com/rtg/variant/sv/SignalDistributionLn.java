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

import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;


/**
 */
public class SignalDistributionLn extends IntegralAbstract implements Signal {

  private final SamCounts mCounts;

  private final Distribution mDistribution;

  private final Distribution mDistributionLn;

  private final int mWindowLo;

  private final int mWindowHi;

  private final String mColumnName;

  /**
   * @param counts SamArray used to construct the signal.
   * @param distr underlying distribution of rates.
   * @param column the name of the column for output.
   */
  public SignalDistributionLn(SamCounts counts, Distribution distr, String column) {
    mCounts = counts;
    mWindowLo = distr.lo();
    mWindowHi = distr.hi();
    final int diameter = mWindowHi - mWindowLo;
    final double[] ln = new double[diameter];
    for (int i = 0; i < diameter; ++i) {
      ln[i] = -Math.log(distr.get(i + mWindowLo));
    }
    mDistribution = distr;
    mDistributionLn = new DistributionArray(mWindowLo, ln);
    mColumnName = column;
    assert globalIntegrity();
  }

  @Override
  public double value(int position) {
    double sum = 0;
    //System.err.println("value position=" + position);
    for (int i = mWindowLo; i < mWindowHi; ++i) {
      final double count = mCounts.count(position, i);
      double s = mDistribution.get(i);
      if (count != 0) {
        s += count * (mDistributionLn.get(i) + Math.log(count) - 1.0);
      }
      if (s < 0.0 && s >= -1.0e-6) {
        s = 0.0; //allow for rounding errors
      }
      //System.err.println("i=" + i + " rho=" + mDistribution.get(i) + " n=" + count + " s=" + s);
      assert s >= 0.0 && Double.isFinite(s) : s + ":" + count + ":" + mDistributionLn.get(i) + ":" + mDistribution.get(i);
      sum += s;
    }
    assert sum >= 0.0 && Double.isFinite(sum);
    return sum;
  }

  @Override
  public void toString(StringBuilder sb) {
    sb.append("DistributionLn:").append(mWindowLo).append(":").append(mWindowHi);
  }

  @Override
  public boolean integrity() {
    Exam.assertNotNull(mCounts);
    Exam.assertNotNull(mDistributionLn);
    Exam.assertTrue(mWindowHi > 0);
    Exam.assertTrue(mWindowLo < 0);
    return true;
  }

  @Override
  public final boolean globalIntegrity() {
    integrity();
    mDistributionLn.globalIntegrity();
    return true;
  }

  @Override
  public String columnLabel() {
    return mColumnName;
  }

}
