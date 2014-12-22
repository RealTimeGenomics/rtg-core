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
package com.rtg.variant.bayes;

import com.rtg.variant.GenomePriorParams;
import com.rtg.variant.bayes.snp.HypothesesSnp;
import com.rtg.variant.util.arithmetic.LogApproximatePossibility;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;

/**
 * Mock Genotype measure for testing
 */
public class MockGenotypeMeasure implements GenotypeMeasure {
  final int mReference;
  final int mBest;
  final double mPosterior;
  final double mNonIdentityPosterior;

  /**
   * Mock measure containing just a bestPosterior
   * @param posterior the posterior to use
   */
  public MockGenotypeMeasure(double posterior) {
    this(0, 0, posterior, 0);
  }

  /**
   *
   * @param reference reference call index
   * @param best best call index
   * @param posterior best call posterior
   * @param nonIdentityPosterior posterior not equal to reference
   */
  public MockGenotypeMeasure(int reference, int best, double posterior, double nonIdentityPosterior) {
    mReference = reference;
    mBest = best;
    mPosterior = posterior;
    mNonIdentityPosterior = nonIdentityPosterior;

  }
  @Override
  public PossibilityArithmetic arithmetic() {
    return LogApproximatePossibility.SINGLETON;
  }

  @Override
  public double measure(int hypothesis) {
    return LogApproximatePossibility.SINGLETON.prob2Poss(1.0 / 10.0);
  }

  @Override
  public int size() {
    throw new UnsupportedOperationException();
  }

  @Override
  public double bestPosterior() {
    return mPosterior;
  }

  @Override
  public double nonIdentityPosterior() {
    return mNonIdentityPosterior;
  }

  @Override
  public int best() {
    return mBest;
  }

  @Override
  public int reference() {
    return mReference;
  }

  @Override
  public Hypotheses<?> hypotheses() {
    return new HypothesesSnp(arithmetic(), GenomePriorParams.builder().create(), false, mReference);
  }
}
