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
package com.rtg.variant.bayes.multisample.forwardbackward;

import java.util.Arrays;

import com.rtg.variant.bayes.AbstractFactor;
import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.Factor;
import com.rtg.variant.bayes.Hypotheses;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;

/**
 * @param <D> description
 */
public final class MutableFactor<D extends Description> extends AbstractFactor<D> {

  private final double[] mValues;

  /**
   * @param hypotheses underlying set of hypotheses.
   * @param arith arithmetic used to represent values.
   * @param values for each hypothesis as probabilities.
   */
  public MutableFactor(final Hypotheses<D> hypotheses, PossibilityArithmetic arith, final double[] values) {
    super(hypotheses, arith);
    mValues = new double[values.length];
    for (int i = 0; i < size(); ++i) {
      mValues[i] = arith.prob2Poss(values[i]);
    }
  }

  /**
   * @param hypotheses underlying set of hypotheses.
   * @param arith arithmetic used to represent values.
   * @param size of hypotheses space
   */
  MutableFactor(Hypotheses<D> hypotheses, PossibilityArithmetic arith, int size) {
    super(hypotheses, arith, size);
    final double zero = arith.zero();
    mValues = new double[size];
    Arrays.fill(mValues, zero);
  }

  /**
   * @param factor factor with which to initialize this factor
   */
  public MutableFactor(final Factor<D> factor) {
    this(factor.hypotheses(), factor.arithmetic(), factor.size());
    for (int i = 0; i < size(); ++i) {
      set(i, factor.p(i));
    }
  }

  @Override
  public double p(final int index) {
    return mValues[index];
  }

  /**
   * @param index selects the hypothesis.
   * @param value what to set the hypothesis value to (a possibility).
   */
  public void set(final int index, final double value) {
    mValues[index] = value;
  }

  @Override
  public boolean isNormalized() {
    return false;
  }
}
