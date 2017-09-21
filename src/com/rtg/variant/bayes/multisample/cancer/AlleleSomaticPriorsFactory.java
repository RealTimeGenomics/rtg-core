/*
 * Copyright (c) 2017. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
 */

package com.rtg.variant.bayes.multisample.cancer;

import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.Hypotheses;
import com.rtg.variant.bayes.HypothesesPowerSet;

/**
 * Factory for producing somatic priors based on given hypotheses.
 */
class AlleleSomaticPriorsFactory<D extends Description> implements SomaticPriorsFactory {

  private final double[][] mQ;
  private final double mMu;

  /**
   * Construct a new factory for the specified hypotheses.
   * @param normalHypotheses to be mutated.
   * @param mu somatic mutation probability
   */
  AlleleSomaticPriorsFactory(final Hypotheses<D> normalHypotheses, final HypothesesPowerSet<D> cancerHypotheses, final double mu) {
    mMu = mu;
    mQ = SomaticPriorsAllele.makeQ(mu, normalHypotheses, cancerHypotheses);
  }

  @Override
  public double[][] somaticQ(final double mu) {
    if (mu != mMu) {
      throw new UnsupportedOperationException("Dynamic mu not supported");
    }
    return mQ;
  }
}
