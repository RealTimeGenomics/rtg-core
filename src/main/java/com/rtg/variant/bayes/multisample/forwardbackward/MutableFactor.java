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
