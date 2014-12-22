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

import com.rtg.util.memo.DoubleFunction;
import com.rtg.util.memo.DoubleMemo;
import com.rtg.util.memo.Function;
import com.rtg.util.memo.Memo;

/**
 */
public final class EntropyRandom {

  private EntropyRandom() { }

  private static final Memo<DoubleMemo> MEMO = new Memo<>(new FunctionEntropy());

  private static class DF implements DoubleFunction {

    private final int mR;

    /**
     * @param r number of alternative counts.
     */
    DF(int r) {
      super();
      mR = r;
    }

    @Override
    public double fn(int i) {
      return entropyLocal(i, mR);
    }

  }

  private static class FunctionEntropy implements Function<DoubleMemo> {
    @Override
    public DoubleMemo fn(int r) {
      return new DoubleMemo(new DF(r));
    }
  }


  /**
   * @param i count for a particular instance.
   * @param r number of alternative counts.
   * @return <code>i * ln(i + 1/r)</code> - part of computing entropy using Laplace estimator.
   */
  static double entropy(final int i, final int r) {
    return MEMO.fn(r).fn(i);
  }

  /**
   * @param i count for a particular instance.
   * @param r number of alternative counts.
   * @return <code>i * ln(i + 1/r)</code> - part of computing entropy using Laplace estimator.
   */
  static double entropyLocal(final int i, final int r) {
    assert i >= 0;
    assert r >= 1;
    return i * Math.log(i + 1.0 / r);
  }

  /**
   * @param total sum of counts.
   * @param counts of individual categories.
   * @param prior the log of the prior.
   * @return return ln of posterior for random distribution.
   */
  static double randomPosteriorLocal(final int total, final int[] counts, double prior) {
    double sum = -entropy(total, 1);
    for (final int c : counts) {
      sum += entropy(c, counts.length);
    }
    assert !Double.isInfinite(sum) && !Double.isNaN(sum);
    return sum + prior;
  }

}
