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
    assert Double.isFinite(sum);
    return sum + prior;
  }

}
