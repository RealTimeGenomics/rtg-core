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

import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;
import com.rtg.variant.bayes.Code;
import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.Hypotheses;

/**
 * Mutate haploid or diploid priors and generate probabilities for each mutation.
 */
abstract class SomaticPriors<D extends Description> extends IntegralAbstract {

  /**
   * @param mu probability of a somatic mutation.
   * @param ref the reference hypothesis.
   * @param priors unnormalized probabilities of transitions.
   * @return normalized transitions including 1-mutation for the reference.
   */
  static double[] mutationNormalize(final double mu, final int ref, final double[] priors) {
    assert Exam.assertDistribution(priors);
    final double[] norm = new double[priors.length];
    double sum = 0.0;
    for (int k = 0; k < priors.length; ++k) {
      if (k != ref) {
        final double d = priors[k] * mu;
        norm[k] = d;
        sum += d;
      }
    }
    norm[ref] = 1.0 - sum;
    return norm;
  }

  /**
   * @param mu probability of a somatic mutation.
   * @param loh probability that there will be a loss of heterozygosity.
   * @param hypotheses the set of current hypotheses.
   * @param initialPriors probabilities of transitions between haploid hypotheses (assumed to be normalized).
   * @return probabilities of somatic transitions between possibly diploid hypotheses.
   * @param <D> description type
   */
  static <D extends Description> double[][] makeQ(final double mu, double loh, final Hypotheses<D> hypotheses, double[][] initialPriors) {
    final int length = hypotheses.size();
    final double[][] q = new double[length][length];
    new SomaticPriors<D>(hypotheses, mu, loh, initialPriors) {
      @Override
      void update(final int k, final int j, final double probability) {
        q[k][j] += probability;
      }
    }.update();
    return q;
  }

  protected final Hypotheses<?> mHypotheses;
  private final double mLoh;
  private final double mMutation;
  private final double[][] mInitialPriors;

  /**
   * @param hypotheses to be mutated.
   * @param mu probability of a single mutation.
   * @param loh probability of loss of heterozygosity.
   * @param initialPriors initial unnormalized haploid transition probabilities.
   */
  SomaticPriors(Hypotheses<D> hypotheses, double mu, double loh, double[][] initialPriors) {
    mHypotheses = hypotheses;
    mLoh = loh;
    mMutation = mu;
    mInitialPriors = initialPriors;
    assert integrity();
    assert globalIntegrity();
  }

  /**
   * Generate all mutations and their probabilities.
   */
  void update() {
    if (mHypotheses.haploid()) {
      updateHaploid();
    } else {
      if (mLoh < 1.0) {
        updateDiploid();
      }
      if (mLoh > 0.0) {
        updateLoh();
      }
    }
  }

  private void updateHaploid() {
    for (int k = 0; k < mHypotheses.size(); ++k) {
      final double[] mut = mutant(k);
      for (int i = 0; i < mHypotheses.description().size(); ++i) {
        update(k, i, mut[i]);
      }
    }
  }

  private void updateDiploid() {
    final int size = mHypotheses.description().size();
    final double[][] mut = new double[size][];
    for (int i = 0; i < size; ++i) {
      mut[i] = mutant(i);
    }
    final double notLoh = 1.0 - mLoh;
    final Code code = mHypotheses.code();
    for (int k = 0; k < mHypotheses.size(); ++k) {
      final int a = code.a(k);
      final int b = code.bc(k);
      final double[] mut0 = mut[a];
      final double[] mut1 = mut[b];
      for (int i = 0; i < mHypotheses.description().size(); ++i) {
        final double prob0 = mut0[i] * notLoh;
        for (int j = 0; j < mHypotheses.description().size(); ++j) {
          final double prob1 = mut1[j];
          final int k2 = code.code(i, j);
          update(k, k2, prob0 * prob1);
        }
      }
    }
  }

  void updateLoh() {
    final double loh = mLoh * 0.5;
    final int size = mHypotheses.description().size();
    final double[][] mut = new double[size][];
    for (int i = 0; i < size; ++i) {
      mut[i] = mutant(i);
    }
    final Code code = mHypotheses.code();
    for (int k = 0; k < mHypotheses.size(); ++k) {
      final int c0 = code.a(k);
      final double[] mut0 = mut[c0];
      for (int j = 0; j < mHypotheses.description().size(); ++j) {
        final double prob = mut0[j];
        update(k, code.code(j, j), loh * prob);
      }

      final int c2 = code.bc(k);
      final double[] mut2 = mut[c2];
      for (int j = 0; j < mHypotheses.description().size(); ++j) {
        final double prob = mut2[j];
        update(k, code.code(j, j), loh * prob);
      }
    }
  }

  /**
   * Called for each mutation.
   * @param key1 the original name of category.
   * @param key2 the mutated name of category.
   * @param probability of the mutation.
   */
  abstract void update(final int key1, final int key2, final double probability);

  /**
   * Compute the transition probability from k to each of the allowed hypotheses.
   * @param k the original hypothesis
   * @return the transition probabilities.
   */
  double[] mutant(int k) {
    return mutationNormalize(mMutation, k, mInitialPriors[k]);
  }

  @Override
  public final boolean globalIntegrity() {
    integrity();
    final int size = mHypotheses.description().size();
    for (int i = 0; i < size; ++i) {
      final double[] pr = mInitialPriors[i];
      Exam.assertEquals(size, pr.length);
      for (int j = 0; j < size; ++j) {
        final double pv = pr[j];
        Exam.assertTrue(0.0 <= pv && pv <= 1.0);
      }
    }
    return true;
  }

  @Override
  public final boolean integrity() {
    Exam.assertTrue(0.0 <= mMutation && mMutation <= 1.0);
    final int size = mHypotheses.description().size();
    Exam.assertEquals(size, mInitialPriors.length);
    return true;
  }
}
