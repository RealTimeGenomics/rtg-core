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

  private final Hypotheses<D> mNormalHypotheses;

  /**
   * Construct a new factory for the specified hypotheses.
   * @param normalHypotheses normal hypotheses
   */
  AlleleSomaticPriorsFactory(final Hypotheses<D> normalHypotheses) {
    mNormalHypotheses = normalHypotheses;
  }

  @Override
  public double[][] somaticQ(final double mu) {
    return SomaticPriorsAllele.makeQ(mu, mNormalHypotheses, new HypothesesPowerSet<>(mNormalHypotheses.description(), mNormalHypotheses.arithmetic(), mNormalHypotheses.reference()));
  }
}
