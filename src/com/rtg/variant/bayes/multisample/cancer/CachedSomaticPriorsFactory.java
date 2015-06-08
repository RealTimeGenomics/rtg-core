/*
 * Copyright (c) 2015. Real Time Genomics Limited.
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

import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.Hypotheses;

/**
 * Factory for producing somatic priors based on given hypotheses.
 */
class CachedSomaticPriorsFactory<D extends Description> extends SomaticPriorsFactory<D> {

  private final double[][][] mCache = new double[-Double.MIN_EXPONENT + 1][][];

  /**
   * Construct a new factory for the specified hypotheses.
   * @param hypotheses to be mutated.
   * @param loh probability of loss of heterozygosity.
   * @param initialPriors initial unnormalized haploid transition probabilities.
   */
  CachedSomaticPriorsFactory(final Hypotheses<D> hypotheses, final double loh, final double[][] initialPriors) {
    super(hypotheses, loh, initialPriors);
  }

  /**
   * Construct a new factory for the specified hypotheses.
   * @param hypotheses to be mutated.
   * @param loh probability of loss of heterozygosity.
   */
  CachedSomaticPriorsFactory(final Hypotheses<D> hypotheses, final double loh) {
    super(hypotheses, loh);
  }

  /**
   * Compute the somatic Q prior matrix for the specified mutation rate.
   * @param mutation probability of a somatic mutation.
   * @return probabilities of somatic transitions between possibly diploid hypotheses.
   */
  double[][] somaticQ(final double mutation) {
    // Binning based on powers of 2, could be precomputed to avoid synchronized
    final int exponent = Math.getExponent(mutation);
    assert exponent <= 1;
    final int index = 1 - exponent;
    assert index >= 0 && index <= mCache.length;
    synchronized (mCache) {
      if (mCache[index] == null) {
        final double binnedMutation = Math.scalb(1, exponent);
        assert binnedMutation <= mutation && 2 * binnedMutation >= mutation;
        mCache[index] = super.somaticQ(binnedMutation);
      }
      return mCache[index];
    }
  }


}
