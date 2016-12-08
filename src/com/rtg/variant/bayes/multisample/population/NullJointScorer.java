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

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.ModelInterface;
import com.rtg.variant.bayes.multisample.HypothesisScore;
import com.rtg.variant.bayes.multisample.HypothesisScores;
import com.rtg.variant.bayes.multisample.MultisampleJointScorer;
import com.rtg.variant.bayes.multisample.PriorContainer;
import com.rtg.variant.bayes.multisample.forwardbackward.BContainer;
import com.rtg.variant.bayes.snp.HypothesesPrior;

/**
 * Records results of one iteration of EM algorithm for population estimators.
 * Makes same estimate as started with (useful for testing).
 */
@TestClass("com.rtg.variant.bayes.multisample.population.NullEstimatorTest")
class NullJointScorer implements MultisampleJointScorer {

  NullJointScorer() { }

  @Override
  public <D extends Description, T extends HypothesesPrior<D>> HypothesisScores getBestScores(List<ModelInterface<?>> models, PriorContainer<T> priorContainer) {
    final HypothesisScore[] calls = new HypothesisScore[models.size()];
    for (int i = 0; i < models.size(); ++i) {
      final ModelInterface<?> model = models.get(i);
      final HypothesisScore call = model.best(priorContainer.getHypotheses().get(model));
      calls[i] = call;
    }
    return new HypothesisScores(calls, false, 0, null);
  }

  @Override
  public BContainer[] makeInitialBs(List<ModelInterface<?>> models) {
    return null;
  }
}
