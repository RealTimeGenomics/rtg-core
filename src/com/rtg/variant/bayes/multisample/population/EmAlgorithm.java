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
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.SpyHistogram;
import com.rtg.util.format.FormatInteger;
import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.ModelInterface;
import com.rtg.variant.bayes.multisample.HypothesisScore;
import com.rtg.variant.bayes.multisample.HypothesisScores;
import com.rtg.variant.bayes.multisample.MultisampleJointScorer;
import com.rtg.variant.bayes.multisample.PriorContainer;
import com.rtg.variant.bayes.multisample.forwardbackward.BContainer;
import com.rtg.variant.bayes.snp.HypothesesPrior;

/**
 * Re-estimates priors using simple Laplace estimator on all hypotheses.
 * Assumes mixture of haploid and diploid models.
 */
@TestClass({"com.rtg.variant.bayes.multisample.population.HwEstimatorTest", "com.rtg.variant.bayes.multisample.population.EmAlgorithmTest"})
public class EmAlgorithm implements MultisampleJointScorer {

  private static final int ITERATION_HIST_SIZE = 50;
  private static final SpyHistogram ITERATION_HIST = new SpyHistogram("EmAlgorithm iterations", ITERATION_HIST_SIZE);

  private static final FormatInteger INT_FORMAT_1 = new FormatInteger(1);
  private static final FormatInteger INT_FORMAT_3 = new FormatInteger(3);
  private static final int ARBITRARY_LOG_TRIGGER = 50;

  /** Default maximum number of iterations */
  public static final int DEFAULT_MAX_ITERATIONS = 50;

  static String callToString(final int size, final HypothesisScore[] calls) {
    final FormatInteger format;
    final String separator;
    if (size <= 10) {
      format = INT_FORMAT_1;
      separator = " ";
    } else {
      format = INT_FORMAT_3;
      separator = " | ";
    }
    final StringBuilder sb = new StringBuilder();
    for (int i = 0; i < calls.length; i++) {
      final int call = calls[i].hypothesis();
      sb.append(format.format(call));
      if (i % 10 == 9) {
        sb.append(separator);
      }
    }
    return sb.toString();
  }

  static synchronized void updateIterationHist(int itcount) {
    ITERATION_HIST.increment(itcount);
  }

  private final Estimator mEstimator;

  private final int mMaxIterations;

  /**
   * @param estimator used for estimating the new priors.
   * @param maxIterations the maximum number of iterations that EM will attempt
   */
  public EmAlgorithm(Estimator estimator, int maxIterations) {
    mEstimator = estimator;
    mMaxIterations = maxIterations < 0 ? DEFAULT_MAX_ITERATIONS : maxIterations;
    //System.err.println("emAlgo: " + mMaxIterations);
  }


  @Override
  public <D extends Description, T extends HypothesesPrior<D>> HypothesisScores getBestScores(List<ModelInterface<?>> models, PriorContainer<T> priorContainer) {
    int iterations = 0;
    EmResult<HypothesesPrior<D>> last = mEstimator.estimate(models, priorContainer);
    //System.err.println(callToString(mInitialDiploid.size(), last.calls()));
    //System.err.println(priorsToString(last.haploid()));
    //System.err.println(priorsToString(last.diploid()));
    //int count = 0;
    while (iterations < mMaxIterations) {
      final EmResult<HypothesesPrior<D>> next = mEstimator.estimate(models, last.getPriorContainer());
      final int difference = last.difference(next.calls());
      //System.err.println("diff=" + difference);
      //System.err.println(callToString(mInitialDiploid.size(), next.calls()));
      //System.err.println(priorsToString(next.haploid()));
      //System.err.println(priorsToString(next.diploid()));
      iterations++;
      last = next;
      if (difference == 0) {
        break;
      }
    }
    updateIterationHist(iterations);
    if (iterations == mMaxIterations && mMaxIterations >= ARBITRARY_LOG_TRIGGER) {
      Diagnostic.userLog("EmAlgorithm exceeded " + mMaxIterations + " iterations");
      Diagnostic.developerLog(ITERATION_HIST.toString());
    }
    return last.calls();
  }

  @Override
  public BContainer[] makeInitialBs(List<ModelInterface<?>> models) {
    return null;
  }
}
