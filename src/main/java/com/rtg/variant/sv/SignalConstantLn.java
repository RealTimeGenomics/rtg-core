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
 * Like SignalDistributionLn but optimized to operate on DistributionConstant or DistributionStep distributions.
 */
public class SignalConstantLn extends IntegralAbstract implements Signal {

  private final SamCounts mCounts;

  private final double mLn1;
  private final double mRate1;

  private final int mOffset;

  private final double mLn2;
  private final double mRate2;

  private final int mWindowLo;

  private final int mWindowHi;


  private final String mColumnName;

  /**
   * @param counts SamArray used to construct the signal.
   * @param distr underlying distribution of rates.
   * @param column the name of the column for output.
   */
  SignalConstantLn(SamCounts counts, DistributionConstant distr, String column) {
    mCounts = counts;
    mWindowLo = distr.lo();
    mWindowHi = distr.hi();
    mRate1 = distr.getConstant();
    mLn1 = -Math.log(mRate1);
    mOffset = mWindowHi;
    mRate2 = 0;
    mLn2 = 0;
    mColumnName = column;
    assert globalIntegrity();
  }

  /**
   * @param counts SamArray used to construct the signal.
   * @param distr underlying distribution of rates.
   * @param column the name of the column for output.
   */
  SignalConstantLn(SamCounts counts, DistributionStep distr, String column) {
    mCounts = counts;
    mWindowLo = distr.lo();
    mWindowHi = distr.hi();
    mRate1 = distr.getRate1();
    mLn1 = -Math.log(mRate1);
    mOffset = distr.getOffset() + 1;
    mRate2 = distr.getRate2();
    mLn2 = -Math.log(mRate2);
    mColumnName = column;
    assert globalIntegrity();
  }

  @Override
  public double value(int position) {
    final int diameter1 = mOffset - mWindowLo;
    final double countTot1 = mCounts.count(position, mWindowLo, mOffset);
    double s = diameter1 * mRate1;
    if (countTot1 > 0) {
      s += countTot1 * (mLn1 - 1.0) + mCounts.sumLn(position, mWindowLo, mOffset);
    }

    final int diameter2 = mWindowHi - mOffset;
    if (diameter2 > 0) {
      final double countTot2 = mCounts.count(position, mOffset, mWindowHi);
      s += diameter2 * mRate2;
      if (countTot2 > 0) {
        s += countTot2 * (mLn2 - 1.0) + mCounts.sumLn(position, mOffset, mWindowHi);
      }
    }
    if (s < 0.0 && s >= -1.0e-6) {
      s = 0.0; //allow for rounding errors
    }
    return s;
  }

  @Override
  public void toString(StringBuilder sb) {
    sb.append("DistributionLn:").append(mWindowLo).append(":").append(mWindowHi);
  }

  @Override
  public boolean integrity() {
    Exam.assertNotNull(mCounts);
    Exam.assertTrue(mWindowLo < 0);
    Exam.assertTrue(mWindowHi > 0);
    return true;
  }

  @Override
  public final boolean globalIntegrity() {
    integrity();
    return true;
  }

  @Override
  public String columnLabel() {
    return mColumnName;
  }

}
