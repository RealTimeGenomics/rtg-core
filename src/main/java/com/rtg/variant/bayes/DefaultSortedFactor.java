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

import java.util.Arrays;
import java.util.Comparator;

import com.rtg.variant.util.arithmetic.PossibilityArithmetic;

/**
 * A sorted factor backed by another factor.
 * @param <D> description backing the factor (the scope).
 */
public class DefaultSortedFactor<D extends Description> implements SortedFactor<D> {

  private final Factor<D> mFactor;
  private final Integer[] mRankToCode;

  private static final class FactorComparator implements Comparator<Integer> {

    private final Factor<?> mFactor;

    private FactorComparator(final Factor<?> factor) {
      mFactor = factor;
    }

    @Override
    public int compare(final Integer a, final Integer b) {
      final double p = mFactor.p(a);
      final double q = mFactor.p(b);
      if (mFactor.arithmetic().gt(p, q)) {
        return -1;
      } else if (mFactor.arithmetic().gt(q, p)) {
        return 1;
      } else {
        return 0;
      }
    }
  }

  private static Integer[] sort(final Factor<?> factor) {
    // Set up identity map
    final Integer[] rankToCode = new Integer[factor.size()];
    for (int k = 0; k < rankToCode.length; ++k) {
      rankToCode[k] = k;
    }
    // Sort with respect to factor
    Arrays.sort(rankToCode, new FactorComparator(factor));
    return rankToCode;
  }

  DefaultSortedFactor(final Factor<D> factor) {
    mFactor = factor;
    mRankToCode = sort(factor);
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
    return mFactor.p(code);
  }

  @Override
  public boolean isNormalized() {
    return mFactor.isNormalized();
  }

  @Override
  public Factor<D> normalize() {
    throw new UnsupportedOperationException("Normalize before sorting");
  }

  @Override
  public int hypothesis(int rank) {
    return mRankToCode[rank];
  }

  @Override
  public double value(int rank) {
    return mFactor.p(hypothesis(rank));
  }
}
