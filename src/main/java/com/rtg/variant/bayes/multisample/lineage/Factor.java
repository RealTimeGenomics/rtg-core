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
