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
  public boolean globalIntegrity() {
    integrity();
    return true;
  }

  @Override
  public String columnLabel() {
    return mColumnName;
  }

}
