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
import com.rtg.variant.bayes.multisample.PriorContainer;
import com.rtg.variant.bayes.snp.HypothesesPrior;

/**
 */
public interface Estimator {

  /**
   * @param <D> type of description
   * @param <T> type of hypotheses
   * @param models which can be interrogated for calls.
   * @param priorContainer box of prior related things
   * @return the result of one iteration of the EM algorithm.
   */
  <D extends Description, T extends HypothesesPrior<D>> EmResult<HypothesesPrior<D>> estimate(final List<ModelInterface<?>> models, PriorContainer<T> priorContainer);
}
