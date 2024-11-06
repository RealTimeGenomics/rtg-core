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
package com.rtg.variant.bayes.multisample.lineage;

import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;

/**
 * Common parts of factor implementation.
 */
@TestClass("com.rtg.variant.bayes.multisample.lineage.DefaultFactorTest")
public abstract class AbstractFactor implements Factor {

  private final PossibilityArithmetic mArithmetic;

  /**
   * Factor.
   * @param arithmetic the arithmetic
   */
  public AbstractFactor(final PossibilityArithmetic arithmetic) {
    mArithmetic = arithmetic;
  }

  @Override
  public PossibilityArithmetic arithmetic() {
    return mArithmetic;
  }

  protected void checkScope(final Map<Variable, Integer> values) {
    if (!values.keySet().equals(scope())) {
      throw new IllegalArgumentException("Scope error");
    }
  }

  @Override
  public abstract double p(final Map<Variable, Integer> values);

  @Override
  public abstract Set<Variable> scope();

  @Override
  public Factor multiply(final Factor other) {
    return DefaultFactor.asDefault(this).multiply(other);
  }

  @Override
  public Factor sumOut(final Variable variable) {
    final HashSet<Variable> keep = new HashSet<>(scope());
    if (!keep.remove(variable)) {
      throw new IllegalArgumentException(variable + " not in scope");
    }
    return marginal(keep);
  }

  @Override
  public Factor marginal(final Set<Variable> variablesToKeep) {
    return DefaultFactor.asDefault(this).marginal(variablesToKeep);
  }

  @Override
  public Factor condition(final Map<Variable, Integer> conditions) {
    return DefaultFactor.asDefault(this).condition(conditions);
  }
}
