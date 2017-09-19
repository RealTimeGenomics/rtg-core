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

import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.Hypotheses;

/**
 * Factory for producing somatic priors based on given hypotheses.  This implementation
 * only computes the Q matrix for powers of 2 of the somatic mutation rate.  Thus, given
 * mu next lower power of 2 is used for the actual calculation.  These powers of 2 values
 * are remembered so they can be reused.
 */
class CachedSomaticPriorsFactory<D extends Description> extends SomaticPriorsFactory<D> {

  // The +2 below is needed to handle mu < Double.MIN_NORMAL (i.e. 0 and subnormal numbers).
  // Cf. Math.exponent()
  private final double[][][] mCache = new double[-Double.MIN_EXPONENT + 2][][];

  /**
   * Construct a new factory for the specified hypotheses.
   * @param hypotheses to be mutated.
   * @param loh probability of loss of heterozygosity.
   */
  CachedSomaticPriorsFactory(final Hypotheses<D> hypotheses, final double loh) {
    super(hypotheses, loh);
  }

  @Override
  double[][] somaticQ(final double mu) {
    // Binning based on powers of 2, could be precomputed to avoid synchronized
    final int exponent = Math.getExponent(mu);
    assert exponent <= 1;
    final int index = -exponent;
    assert index >= 0 && index < mCache.length : String.valueOf(index);
    synchronized (mCache) {
      if (mCache[index] == null) {
        final double binnedMutation = Math.scalb(1, exponent);
        assert binnedMutation <= mu && binnedMutation >= 0.5 * mu;
        mCache[index] = super.somaticQ(binnedMutation);
      }
      return mCache[index];
    }
  }
}
