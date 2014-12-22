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
package com.rtg.variant.bayes.complex;

import com.rtg.util.integrity.Exam;
import com.rtg.util.memo.DoubleFunction;
import com.rtg.util.memo.DoubleMemo;

/**
 */
public class KappaMemo extends AbstractKappa {

  private final AbstractKappa mKappa;

  private final DoubleMemo mPiMemo;

  private final DoubleMemo mPiSumMemo;

  /**
   * @param kappa underlying implementation being memoized.
   */
  public KappaMemo(AbstractKappa kappa) {
    mKappa = kappa;
    mPiMemo = new DoubleMemo(new DoubleFunction() {
      @Override
      public double fn(int i) {
        return mKappa.pi(i);
      }
    });
    mPiSumMemo = new DoubleMemo(new DoubleFunction() {
      @Override
      public double fn(int i) {
        return mKappa.piSum(i);
      }
    });
  }

  @Override
  double pi(int i) {
    return mPiMemo.fn(i);
  }

  @Override
  double piSum(int i) {
    return mPiSumMemo.fn(i);
  }

  @Override
  public void toString(StringBuilder sb) {
    sb.append("KappaMemo ");
    mKappa.toString(sb);
  }

  @Override
  public boolean globalIntegrity() {
    integrity();
    mKappa.globalIntegrity();
    mPiMemo.globalIntegrity();
    mPiSumMemo.globalIntegrity();
    return true;
  }

  @Override
  public boolean integrity() {
    Exam.assertNotNull(mKappa);
    Exam.assertNotNull(mPiMemo);
    Exam.assertNotNull(mPiSumMemo);
    return true;
  }

}
