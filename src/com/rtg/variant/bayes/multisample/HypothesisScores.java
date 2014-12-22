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

import com.rtg.variant.bayes.multisample.forwardbackward.BContainer;

/**
 * Wraps up some joint caller scoring
 */
public class HypothesisScores {

  final HypothesisScore[] mScores;
  final boolean mInteresting;
  final double mNonIdentityPosterior;
  final BContainer[] mBs;

  /**
   * Wrapper for a set of raw multi-sample calls
   * @param scores the best scoring hypotheses for each model
   * @param interesting true if the overall call is interesting
   * @param nonIdentityPosterior the overall non identity posterior
   * @param bs new value for <code>B's</code> to be passed back to EM algorithm
   */
  public HypothesisScores(HypothesisScore[] scores, boolean interesting, double nonIdentityPosterior, BContainer[] bs) {
    mScores = scores;
    mInteresting = interesting;
    mNonIdentityPosterior = nonIdentityPosterior;
    mBs = bs;
  }

  public HypothesisScore[] getScores() {
    return mScores;
  }

  public BContainer[] getBs() {
    return mBs;
  }

  public boolean isInteresting() {
    return mInteresting;
  }

  public double getNonIdentityPosterior() {
    return mNonIdentityPosterior;
  }

}
