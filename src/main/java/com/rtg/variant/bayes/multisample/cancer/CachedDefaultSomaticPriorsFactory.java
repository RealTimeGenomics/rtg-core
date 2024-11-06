/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
class CachedDefaultSomaticPriorsFactory<D extends Description> extends DefaultSomaticPriorsFactory<D> {

  // The +2 below is needed to handle mu < Double.MIN_NORMAL (i.e. 0 and subnormal numbers).
  // Cf. Math.exponent()
  private final double[][][] mCache = new double[-Double.MIN_EXPONENT + 2][][];

  /**
   * Construct a new factory for the specified hypotheses.
   * @param hypotheses to be mutated.
   * @param loh probability of loss of heterozygosity.
   */
  CachedDefaultSomaticPriorsFactory(final Hypotheses<D> hypotheses, final double loh) {
    super(hypotheses, loh);
  }

  @Override
  public double[][] somaticQ(final double mu) {
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
