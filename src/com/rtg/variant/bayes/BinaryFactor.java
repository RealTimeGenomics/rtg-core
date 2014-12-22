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

import com.rtg.variant.util.arithmetic.PossibilityArithmetic;

/**
 * Factor for a binary hypothesis.
 * @param <D> description
 */
public class BinaryFactor<D extends Description> extends AbstractFactor<D> {

  private final double[] mPoss;

  /**
   * @param hypotheses underlying hypotheses.
   * @param arith arithmetic used to represent values.
   * @param prob probability to associate with hypothesis 1 (true).
   */
  public BinaryFactor(final Hypotheses<D> hypotheses, final PossibilityArithmetic arith, final double prob) {
    super(hypotheses, arith);
    if (prob < 0 || prob > 1 || hypotheses.size() != 2) {
      throw new IllegalArgumentException();
    }
    mPoss = new double[] {arithmetic().prob2Poss(1 - prob), arithmetic().prob2Poss(prob)};
  }

  @Override
  public double p(final int index) {
    return mPoss[index];
  }

  @Override
  public boolean isNormalized() {
    return true;
  }
}
