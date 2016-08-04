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

package com.rtg.variant.realign;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;
import com.rtg.variant.util.arithmetic.LogApproximatePossibility;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;
import com.rtg.variant.util.arithmetic.SimplePossibility;

/**
 * Common code for scoring that first tries with simple doubles and if this
 * underflows tries again using logs (which are slower).
 */
@TestClass(value = "com.rtg.variant.realign.ScoreFastUnderflowTest")
public abstract class AbstractAllPathsFastUnderflow extends IntegralAbstract implements AllPaths {

  private final RealignParams mParams;

  private AllPaths mFast = null;
  private AllPaths mSlowSure = null;
  private AllPaths mUseThis = null;

  /**
   * @param params the machine error model and related parameters.
   */
  public AbstractAllPathsFastUnderflow(final RealignParams params) {
    mParams = params;
  }

  /**
   * Make a new scoring matrix.
   * @param arith helper object that does the arithmetic so that this code can be independent of the representation.
   * @param params the gap probabilities to use.
   * @return the new matrix
   */
  protected abstract AllPaths makeMatrix(final PossibilityArithmetic arith, final RealignParams params);

  @Override
  public void setEnv(Environment env) {
    if (mFast == null) {
      mFast = makeMatrix(SimplePossibility.SINGLETON, mParams);
    }
    mFast.setEnv(env);
    if (mFast.underflow()) {
      //System.err.println("Underflow");
      if (mSlowSure == null) {
        mSlowSure = makeMatrix(LogApproximatePossibility.SINGLETON, mParams);
      }
      mSlowSure.setEnv(env);
      mUseThis = mSlowSure;
    } else {
      mUseThis = mFast;
    }
  }

  @Override
  public PossibilityArithmetic arithmetic() {
    return mUseThis.arithmetic();
  }

  @Override
  public double totalScoreLn() {
    return mUseThis.totalScoreLn();
  }

  @Override
  public double totalScore() {
    return mUseThis.totalScore();
  }

  @Override
  public double total() {
    return mUseThis.total();
  }

  @Override
  public boolean underflow() {
    return false;
  }

  @Override
  public String toString() {
    return mUseThis.toString();
  }

  @Override
  public boolean integrity() {
    Exam.assertNotNull(mParams);
    Exam.assertNotNull(mFast);
    if (mUseThis == null) {
      Exam.assertTrue(mSlowSure == null);
    } else if (mUseThis == mFast) {
      Exam.assertTrue(!mFast.underflow());
    } else {
      Exam.assertTrue(mFast.underflow());
      Exam.assertNotNull(mSlowSure);
      Exam.assertTrue(mUseThis == mSlowSure);
    }
    return true;
  }
}
