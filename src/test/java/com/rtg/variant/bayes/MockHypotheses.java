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

import com.rtg.variant.bayes.snp.HypothesesPrior;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;

/**
 * @param <D> description type
 */
public class MockHypotheses<D extends Description> extends HypothesesPrior<D> {
  private final double[] mPriors;

  /**
   * @param description of the haploid hypotheses.
   * @param arithmetic used for calculations.
   * @param haploid true iff to construct a haploid set of hypotheses.
   * @param priors prior probabilities.
   * @param reference reference hypothesis
   */
  public MockHypotheses(D description, PossibilityArithmetic arithmetic, boolean haploid, double[] priors, int reference) {
    super(description, arithmetic, haploid, reference);
    if (priors != null) {
      mPriors = new double[size()];
      for (int i = 0; i < size(); ++i) {
        mPriors[i] = arithmetic.prob2Poss(priors[i]);
      }
    } else {
      mPriors = null;
    }
  }

  @Override
  public double p(int hyp) {
    return mPriors[hyp];
  }

}
