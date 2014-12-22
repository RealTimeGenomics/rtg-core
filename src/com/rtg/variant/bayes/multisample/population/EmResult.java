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

import com.rtg.variant.bayes.multisample.HaploidDiploidHypotheses;
import com.rtg.variant.bayes.multisample.HypothesisScore;
import com.rtg.variant.bayes.multisample.HypothesisScores;
import com.rtg.variant.bayes.multisample.PriorContainer;
import com.rtg.variant.bayes.multisample.forwardbackward.BContainer;
import com.rtg.variant.bayes.snp.HypothesesPrior;

/**
 * Result of one iteration of EM algorithm.
 * @param <T> type of hypotheses
 */
public class EmResult<T extends HypothesesPrior<?>> {

  private final HypothesisScores mCalls;

  private final PriorContainer<T> mPriorContainer;


  /**
   * @param hypotheses the haploid diploid hypotheses
   * @param calls made by the models.
   * @param bs B values for all samples
   */
  public EmResult(HaploidDiploidHypotheses<T> hypotheses, HypothesisScores calls, BContainer[] bs) {
    mPriorContainer = new PriorContainer<>(hypotheses, bs);
    mCalls = calls;
  }

  PriorContainer<T> getPriorContainer() {
    return mPriorContainer;
  }

  int difference(final HypothesisScores calls) {
    return difference(calls.getScores(), mCalls.getScores());
  }

  static int difference(final HypothesisScore[] calls0, final HypothesisScore[] calls1) {
    assert calls0.length == calls1.length;
    int count = 0;
    for (int i = 0; i < calls0.length; i++) {
      if (calls0[i] == null) { // E.g. female on Y chromosome
        if (calls1[i] != null) { // Can't suddenly change across em iterations.
          throw new RuntimeException();
        }
      } else if (calls0[i].hypothesis() != calls1[i].hypothesis()) {
        count++;
      }
    }
    return count;
  }

  /**
   * Get calls.
   * @return Returns the calls.
   */
  HypothesisScores calls() {
    return mCalls;
  }

}
