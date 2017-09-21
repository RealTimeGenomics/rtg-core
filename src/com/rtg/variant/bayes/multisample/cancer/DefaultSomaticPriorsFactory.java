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

package com.rtg.variant.bayes.multisample.cancer;

import java.util.Arrays;

import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.Hypotheses;

/**
 * Factory for producing somatic priors based on given hypotheses.
 */
class DefaultSomaticPriorsFactory<D extends Description> implements SomaticPriorsFactory {

  private static double[][] defaultUniformPriors(final int size) {
    // Each row is normalized with zero probability for the identity
    final double uniform = 1.0 / (size - 1);
    final double[][] initialPriors = new double[size][size];
    for (int k = 0; k < size; ++k) {
      Arrays.fill(initialPriors[k], uniform);
      initialPriors[k][k] = 0;
    }
    return initialPriors;
  }

  protected final Hypotheses<D> mHypotheses;
  private final double mLoh;
  private final double[][] mInitialPriors;

  /**
   * Construct a new factory for the specified hypotheses.
   * @param hypotheses to be mutated.
   * @param loh probability of loss of heterozygosity.
   * @param initialPriors initial unnormalized haploid transition probabilities.
   */
  DefaultSomaticPriorsFactory(final Hypotheses<D> hypotheses, final double loh, final double[][] initialPriors) {
    mHypotheses = hypotheses;
    mLoh = loh;
    mInitialPriors = initialPriors;
  }

  /**
   * Construct a new factory for the specified hypotheses.
   * @param hypotheses to be mutated.
   * @param loh probability of loss of heterozygosity.
   */
  DefaultSomaticPriorsFactory(final Hypotheses<D> hypotheses, final double loh) {
    this(hypotheses, loh, defaultUniformPriors(hypotheses.description().size()));
  }

  @Override
  public double[][] somaticQ(final double mu) {
    final int length = mHypotheses.size();
    final double[][] q = new double[length][length];
    new SomaticPriors<D>(mHypotheses, mu, mLoh, mInitialPriors) {
      @Override
      void update(final int k, final int j, final double probability) {
        q[k][j] += probability;
      }
    }.update();
    return q;
  }
}
