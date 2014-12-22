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
