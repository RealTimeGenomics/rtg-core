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
