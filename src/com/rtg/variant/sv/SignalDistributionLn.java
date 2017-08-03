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
  public boolean globalIntegrity() {
    integrity();
    mDistributionLn.globalIntegrity();
    return true;
  }

  @Override
  public String columnLabel() {
    return mColumnName;
  }

}
