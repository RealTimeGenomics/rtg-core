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

import com.rtg.util.integrity.Exam;
import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.Hypotheses;

/**
 * @param <D> description type
 */
public abstract class CombinedPriorsComplex<D extends Description> extends CombinedPriors<D> {

  /**
   * @param mutation probability of a somatic mutation.
   * @param ref the reference hypothesis.
   * @param priors unnormalized probabilities of transitions.
   * @return normalized transitions including 1-mutation for the reference.
   */
  static double[] mutationNormalize(final double mutation, final int ref, final double[] priors) {
    assert Exam.assertDistribution(priors);
    final double[] norm = new double[priors.length];
    double sum = 0.0;
    for (int i = 0; i < priors.length; i++) {
      if (i != ref) {
        final double d = priors[i] * mutation;
        norm[i] = d;
        sum += d;
      }
    }
    norm[ref] = 1.0 - sum;
    return norm;
  }

  static double[][] defaultUniformPriors(final int size) {
    final double snpTransition = 1.0 / (size - 1);
    final double[][] initialPriors = new double[size][size];
    for (int k = 0; k < size; k++) {
      Arrays.fill(initialPriors[k], snpTransition);
      initialPriors[k][k] = 0;
    }
    return initialPriors;
  }

  /**
   * @param mutation probability of a somatic mutation.
   * @param loh probability that there will be a loss of heterozygosity.
   * @param hypotheses the set of current hypotheses.
   * @param initialPriors probabilities of transitions between haploid hypotheses (assumed to be normalized).
   * @return probabilities of somatic transitions between possibly diploid hypotheses.
   * @param <D> description type
   */
  static <D extends Description> double[][] makeQ(final double mutation, double loh, final Hypotheses<D> hypotheses, double[][] initialPriors) {
    final int length = hypotheses.size();
    final double[][] q = new double[length][length];
    new CombinedPriorsComplex<D>(hypotheses, mutation, loh, initialPriors) {
      @Override
      void update(int i1, int i2, double probability) {
        q[i1][i2] += probability;
      }
    }.update();
    return q;
  }

  /**
   * Make the Q matrix for the special case where all transitions are equally likely.
   * @param mutation probability of a somatic mutation.
   * @param loh probability that there will be a loss of heterozygosity.
   * @param hypotheses the set of current hypotheses.
   * @return probabilities of somatic transitions between possibly diploid hypotheses.
   * @param <D> description type
   */
  static <D extends Description> double[][] makeQ(final double mutation, double loh, final Hypotheses<D> hypotheses) {
    return makeQ(mutation, loh, hypotheses, defaultUniformPriors(hypotheses.description().size()));
  }




  private final double mMutation;

  private final double[][] mInitialPriors;

  /**
   * @param hypotheses to be mutated.
   * @param mutation probability of a single mutation.
   * @param loh probability of loss of heterozygosity.
   * @param initialPriors initial unnormalized haploid transition probabilities.
   */
  CombinedPriorsComplex(Hypotheses<D> hypotheses, double mutation, double loh, double[][] initialPriors) {
    super(hypotheses, loh);
    mMutation = mutation;
    mInitialPriors = initialPriors;
    assert integrity();
    assert globalIntegrity();
  }

  @Override
  double[] mutant(int k) {
    return mutationNormalize(mMutation, k, mInitialPriors[k]);
  }

  @Override
  public boolean globalIntegrity() {
    integrity();
    final int size = mHypotheses.description().size();
    for (int i = 0; i < size; i++) {
      final double[] pr = mInitialPriors[i];
      Exam.assertEquals(size, pr.length);
      for (int j = 0; j < size; j++) {
        final double pv = pr[j];
        Exam.assertTrue(0.0 <= pv && pv <= 1.0 && !Double.isNaN(pv));
      }
    }
    return true;
  }

  @Override
  public boolean integrity() {
    Exam.assertTrue(0.0 <= mMutation && mMutation <= 1.0 && !Double.isNaN(mMutation));
    final int size = mHypotheses.description().size();
    Exam.assertEquals(size, mInitialPriors.length);
    return true;
  }

}
