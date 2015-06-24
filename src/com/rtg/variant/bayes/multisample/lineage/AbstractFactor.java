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
