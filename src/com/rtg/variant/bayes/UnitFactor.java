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
 * Unit factor.
 * @param <D> description
 */
public class UnitFactor<D extends Description> extends AbstractFactor<D> {

  private final double mUnit;

  /**
   * @param hypotheses underlying hypotheses.
   * @param arith arithmetic used to represent values.
   * @param size number of items in vector, this should be the number of
   * entries in the hypotheses for normal cases, or 1 for the
   * <code>HypothesisNone</code> case
   */
  public UnitFactor(Hypotheses<D> hypotheses, PossibilityArithmetic arith, int size) {
    super(hypotheses, arith, size);
    mUnit = arithmetic().one();
  }

  /**
   * @param hypotheses underlying hypotheses.
   * @param arith arithmetic used to represent values.
   */
  public UnitFactor(Hypotheses<D> hypotheses, PossibilityArithmetic arith) {
    super(hypotheses, arith);
    mUnit = arithmetic().one();
  }

  @Override
  public double p(int index) {
    if (index >= size()) {
      throw new IllegalArgumentException();
    }
    return mUnit;
  }

  @Override
  public boolean isNormalized() {
    return true;
  }
}
