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

import java.util.Map;
import java.util.Set;

import com.rtg.variant.util.arithmetic.PossibilityArithmetic;

/**
 * A function over a set of random variables (possibly empty) to non-negative
 * real numbers.  Subsumes the notion of a probability distribution and
 * conditional probability distribution.
 *
 */
public interface Factor {

  /**
   * Return the arithmetic associated with the factor.
   * @return the arithmetic
   */
  PossibilityArithmetic arithmetic();

  /**
   * Return the affinity from this factor for the hypothesis corresponding
   * to the supplied values.  The returned value is in the arithmetic of the
   * factor and need not be normalized.
   * @param values a value for each variable in the scope of the factor
   * @return non-negative affinity of the hypothesis
   * @throws IllegalArgumentException if the supplied values do not correspond
   * to the scope of the factor.
   */
  double p(Map<Variable, Integer> values);

  /**
   * The set of variables making up the scope of this factor.
   * @return the scope
   */
  Set<Variable> scope();

  /**
   * Form the factor product of this factor with another factor.   The scope of the
   * resulting factor is the union of the scopes of the multiplicands.
   * @param other factor to multiply this facto by
   * @return factor product
   */
  Factor multiply(Factor other);

  /**
   * Produce a new factor that is this factor with the given random variable
   * has been summed out.
   * @param variable variable to sum out
   * @return factor without specified variable
   * @throws IllegalArgumentException if given variable does not exist in this factor.
   */
  Factor sumOut(final Variable variable);

  /**
   * Produce a new factor that is this factor marginalized down to the given variables.
   * The arithmetic of the marginalized factor is the same as this factor.
   * @param variablesToKeep the scope of the marginalized factor
   * @return marginalized factor
   * @throws IllegalArgumentException if a requested variable does not exist in this factor.
   */
  Factor marginal(final Set<Variable> variablesToKeep);

  /**
   * Reduce this factor by conditioning on a set of variables
   * @param conditions the variables values to condition on
   * @return a new factor conditioned appropriately
   */
  Factor condition(Map<Variable, Integer> conditions);
}
