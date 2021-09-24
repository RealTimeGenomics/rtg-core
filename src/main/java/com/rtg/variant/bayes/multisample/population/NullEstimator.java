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

package com.rtg.variant.bayes.multisample.population;

import java.util.List;

import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.ModelInterface;
import com.rtg.variant.bayes.multisample.HaploidDiploidHypotheses;
import com.rtg.variant.bayes.multisample.HypothesisScores;
import com.rtg.variant.bayes.multisample.PriorContainer;
import com.rtg.variant.bayes.snp.HypothesesPrior;

/**
 * Records results of one iteration of EM algorithm for population estimators.
 * Makes same estimate as started with (useful for testing).
 */
class NullEstimator implements Estimator {

  NullEstimator() { }

  @Override
  public <D extends Description, T extends HypothesesPrior<D>> EmResult<HypothesesPrior<D>> estimate(final List<ModelInterface<?>> models, PriorContainer<T> priorContainer) {
    final HypothesisScores allCalls = new NullJointScorer().getBestScores(models, priorContainer);
    final HaploidDiploidHypotheses<T> hypotheses = priorContainer.getHypotheses();
    final HypothesesPrior<D> newHaploid = new HypothesesPrior<>(hypotheses.haploid().description(), hypotheses.haploid().arithmetic(), hypotheses.haploid(), true, hypotheses.haploid().reference());
    final HypothesesPrior<D> newDiploid = new HypothesesPrior<>(hypotheses.diploid().description(), hypotheses.diploid().arithmetic(), hypotheses.diploid(), false, hypotheses.diploid().reference());
    final HaploidDiploidHypotheses<HypothesesPrior<D>> newHapHypDip = new HaploidDiploidHypotheses<>(hypotheses.none(), newHaploid, newDiploid, false, hypotheses.getDescriptionCounts());
    return new EmResult<>(newHapHypDip, allCalls, priorContainer.getBs());
  }

  @Override
  public String toString() {
    return "NullEstimator";
  }

}
