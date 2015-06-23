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
