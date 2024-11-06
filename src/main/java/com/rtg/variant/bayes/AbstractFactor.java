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

package com.rtg.variant.bayes;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;

/**
 * Default implementation of some factor methods.
 * @param <D> description
 */
@TestClass(value = {"com.rtg.variant.bayes.UnitFactorTest"})
public abstract class AbstractFactor<D extends Description> implements Factor<D> {

  protected final Hypotheses<D> mHypotheses;
  protected final PossibilityArithmetic mArithmetic;
  protected final int mSize;

  /**
   * @param hypotheses underlying hypotheses.
   * @param arithmetic arithmetic used to represent values.
   * @param size number of items in vector, this should be the number of entries in the hypotheses for normal cases, or 1 for the <code>HypothesisNone</code> case
   */
  public AbstractFactor(Hypotheses<D> hypotheses, PossibilityArithmetic arithmetic, int size) {
    mHypotheses = hypotheses;
    mArithmetic = arithmetic;
    mSize = size;
  }

  /**
   * @param hypotheses underlying hypotheses.
   * @param arithmetic arithmetic used to represent values.
   */
  public AbstractFactor(Hypotheses<D> hypotheses, PossibilityArithmetic arithmetic) {
    this(hypotheses, arithmetic, hypotheses.size());
  }

  @Override
  public abstract double p(int index);

  @Override
  public Hypotheses<D> hypotheses() {
    return mHypotheses;
  }

  @Override
  public int size() {
    return mSize;
  }

  @Override
  public PossibilityArithmetic arithmetic() {
    return mArithmetic;
  }

  @Override
  public abstract boolean isNormalized();

  @Override
  public Factor<D> normalize() {
    if (isNormalized()) {
      return this;
    }
    return new NormalizedFactor<>(this);
  }

}
