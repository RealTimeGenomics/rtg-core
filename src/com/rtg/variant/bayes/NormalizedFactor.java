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
package com.rtg.variant.bayes;

import com.rtg.variant.util.arithmetic.PossibilityArithmetic;

/**
 * A normalized factor backed by another factor.
 * @param <D> description backing the factor (the scope).
 */
public class NormalizedFactor<D extends Description> implements Factor<D> {

  private final Factor<D> mFactor;
  private final double mInvNormalizer;

  NormalizedFactor(final Factor<D> factor) {
    mFactor = factor;
    final PossibilityArithmetic arithmetic = factor.arithmetic();
    double sum = arithmetic.zero();
    for (int k = 0; k < factor.size(); k++) {
      sum = arithmetic.add(sum, factor.p(k));
    }
    if (arithmetic.isZero(sum)) {
      throw new ArithmeticException();
    }
    mInvNormalizer = arithmetic.prob2Poss(arithmetic.divide(arithmetic.one(), sum));
  }

  @Override
  public Hypotheses<D> hypotheses() {
    return mFactor.hypotheses();
  }

  @Override
  public int size() {
    return mFactor.size();
  }

  @Override
  public PossibilityArithmetic arithmetic() {
    return mFactor.arithmetic();
  }

  @Override
  public double p(int code) {
    return arithmetic().multiply(mFactor.p(code), mInvNormalizer);
  }

  @Override
  public boolean isNormalized() {
    return true;
  }

  @Override
  public Factor<D> normalize() {
    return this;
  }
}
