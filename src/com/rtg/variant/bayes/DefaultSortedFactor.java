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
    for (int k = 0; k < rankToCode.length; k++) {
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
