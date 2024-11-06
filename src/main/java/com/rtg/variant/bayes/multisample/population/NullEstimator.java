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
