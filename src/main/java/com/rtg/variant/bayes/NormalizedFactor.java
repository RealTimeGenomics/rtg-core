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
    for (int k = 0; k < factor.size(); ++k) {
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
