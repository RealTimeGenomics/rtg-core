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

package com.rtg.variant.bayes.multisample;

import java.util.List;

import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.ModelInterface;
import com.rtg.variant.bayes.multisample.forwardbackward.BContainer;
import com.rtg.variant.bayes.snp.HypothesesPrior;


/**
 * Generate best scores for multiple individuals.
 */
public interface MultisampleJointScorer {

  /**
   * Generate best scores for multiple individuals
   *
   * @param <D> the type of the description
   * @param <T> the type of the hypotheses prior
   * @param models input models to call from
   * @param priorContainer container for priors and <code>Bs</code>
   * @return the scores.
   */
  <D extends Description, T extends HypothesesPrior<D>> HypothesisScores getBestScores(List<ModelInterface<?>> models, PriorContainer<T> priorContainer);

  /**
   * If required by scorer create initial B, otherwise null
   * @param models models for all the samples
   * @return initial B values
   */
  BContainer[] makeInitialBs(List<ModelInterface<?>> models);
}
