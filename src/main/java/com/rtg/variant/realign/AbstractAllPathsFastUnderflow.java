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
