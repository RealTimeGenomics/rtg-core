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

import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;
import com.rtg.variant.bayes.Code;
import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.Hypotheses;

/**
 * Mutate haploid or diploid priors and generate probabilities for each mutation.
 */
abstract class CombinedPriors<D extends Description> extends IntegralAbstract {

  protected final Hypotheses<?> mHypotheses;

  private final double mLoh;

  /**
   * @param hypotheses to be mutated.
   * @param loh probability that there will be a loss of heterozygosity.
   */
  CombinedPriors(Hypotheses<?> hypotheses, double loh) {
    mHypotheses = hypotheses;
    mLoh = loh;
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
    for (int k = 0; k < mHypotheses.size(); k++) {
      final double[] mut = mutant(k);
      for (int i = 0; i < mHypotheses.description().size(); i++) {
        update(k, i, mut[i]);
      }
    }
  }

  private void updateDiploid() {
    final int size = mHypotheses.description().size();
    final double[][] mut = new double[size][];
    for (int i = 0; i < size; i++) {
      mut[i] = mutant(i);
    }
    final double notLoh = 1.0 - mLoh;
    final Code code = mHypotheses.code();
    for (int k = 0; k < mHypotheses.size(); k++) {
      final int c0 = code.a(k);
      final int c2 = code.bc(k);
      final double[] mut0 = mut[c0];
      final double[] mut2 = mut[c2];
      for (int i = 0; i < mHypotheses.description().size(); i++) {
        final double prob0 = mut0[i] * notLoh;
        for (int j = 0; j < mHypotheses.description().size(); j++) {
          final double prob2 = mut2[j];
          final int k2 = code.code(i, j);
          // System.out.println("  update[" + k + "][" + k2 + "] = " + prob0 + " * " + prob2 + "  =  " + prob0*prob2);
          update(k, k2, prob0 * prob2);
        }
      }
    }
  }

  void updateLoh() {
    final double loh = mLoh * 0.5;
    final int size = mHypotheses.description().size();
    final double[][] mut = new double[size][];
    for (int i = 0; i < size; i++) {
      mut[i] = mutant(i);
    }
    final Code code = mHypotheses.code();
    for (int k = 0; k < mHypotheses.size(); k++) {
      final int c0 = code.a(k);
      final double[] mut0 = mut[c0];
      for (int j = 0; j < mHypotheses.description().size(); j++) {
        final double prob = mut0[j];
        update(k, code.code(j, j), loh * prob);
      }

      final int c2 = code.bc(k);
      final double[] mut2 = mut[c2];
      for (int j = 0; j < mHypotheses.description().size(); j++) {
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
   * Compute the transition probability from k to each of the allowed codes.
   * @param k the original hypothesis
   * @return the transition probabilities.
   */
  abstract double[] mutant(final int k);

  @Override
  public boolean integrity() {
    Exam.assertNotNull(mHypotheses);
    Exam.assertTrue(0.0 <= mLoh && mLoh <= 1.0 && !Double.isNaN(mLoh));
    return true;
  }


}
