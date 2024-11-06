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
 * A function over a set of random variables (possibly empty) to non-negative
 * real numbers.  Subsumes the notion of a probability distribution and
 * conditional probability distribution.  For an extension to access factor
 * affinities in rank order, see <code>SortedFactor</code>.
 * @param <D> description backing the factor (the scope).
 */
public interface Factor<D extends Description> {

  /**
   * Return the hypotheses corresponding to the factor.  These are
   * the hypotheses from which the integer code used in querying
   * the factor are generated from.
   * @return hypotheses
   */
  Hypotheses<D> hypotheses();

  /**
   * Size of the hypotheses space, except for an empty hypotheses space when
   * the size is 1.  Effectively, one more than the maximum code.
   * @return size of hypotheses space
   */
  int size();

  /**
   * Return the arithmetic associated with the factor.
   * @return the arithmetic
   */
  PossibilityArithmetic arithmetic();

  /**
   * Return the affinity from this factor for the hypothesis corresponding
   * to the supplied code.  The returned value should be a non-negative
   * value, but need not be normalized.  The returned value is in the
   * arithmetic of the factor.
   * @param code code for a hypothesis
   * @return non-negative affinity of hypothesis in factor
   */
  double p(final int code);

  /**
   * Returns true iff this factor is known to be normalized, in the sense that each
   * p-value if a probability and the sum across all hypotheses is 1.
   * This could sometimes return false, even though the factor is normalized; but
   * if it returns true then the factor is definitely normalized.
   * @return true if normalized
   */
  boolean isNormalized();

  /**
   * Return a normalized version of this factor.  That is, the sum of all the
   * affinities in the returned factor will be 1.  If the factor is already
   * normalized, then the same factor may be returned.
   * @return normalized factor
   */
  Factor<D> normalize();
}
