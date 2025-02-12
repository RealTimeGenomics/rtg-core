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

import com.rtg.variant.util.arithmetic.PossibilityArithmetic;

/**
 * Unit factor.
 * @param <D> description
 */
public final class UnitFactor<D extends Description> extends AbstractFactor<D> {

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
