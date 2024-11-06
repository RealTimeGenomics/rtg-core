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

import java.util.Arrays;

import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.Hypotheses;

/**
 * Factory for producing somatic priors based on given hypotheses.
 */
class DefaultSomaticPriorsFactory<D extends Description> implements SomaticPriorsFactory {

  static double[][] defaultUniformPriors(final int size) {
    // Each row is normalized
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
