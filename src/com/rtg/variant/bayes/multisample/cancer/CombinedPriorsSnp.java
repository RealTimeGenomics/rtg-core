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
import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.Hypotheses;

/**
 * @param <D> description type
 */
public abstract class CombinedPriorsSnp<D extends Description> extends CombinedPriors<D> {

  static <D extends Description> double[][] makeQ(final double mutation, double loh, final Hypotheses<D> hypotheses) {
    final int length = hypotheses.size();
    final double[][] q = new double[length][length];
    new CombinedPriorsSnp<D>(hypotheses, mutation, loh) {
      @Override
      void update(int i1, int i2, double probability) {
        q[i1][i2] += probability;
      }
    }.update();
    return q;
  }

  private final double mMutationC;

  private final double mMutationDiv;

  /**
   * @param hypotheses to be mutated.
   * @param mutation probability of a single mutation.
   * @param loh probability that there will be a loss of heterozygosity.
   */
  CombinedPriorsSnp(Hypotheses<?> hypotheses, double mutation, double loh) {
    super(hypotheses, loh);
    mMutationC = 1.0 - mutation;
    mMutationDiv = mutation / (mHypotheses.description().size() - 1);
    assert integrity();
  }

  @Override
  double[] mutant(int k) {
    final double[] mutant = new double[mHypotheses.description().size()];
    for (int i = 0; i < mutant.length; i++) {
      mutant[i] = i == k ? mMutationC : mMutationDiv;
    }
    return mutant;
  }

  @Override
  public boolean integrity() {
    Exam.assertTrue(0.0 <= mMutationC && mMutationC <= 1.0 && !Double.isNaN(mMutationC));
    Exam.assertTrue(0.0 <= mMutationDiv && mMutationDiv <= 1.0 && !Double.isNaN(mMutationDiv));
    return true;
  }

}
