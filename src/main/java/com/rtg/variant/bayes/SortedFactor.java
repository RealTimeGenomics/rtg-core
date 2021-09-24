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

/**
 * A factor allowing retrieval of hypothesis affinity pairs in sorted order.
 * @param <D> description backing the factor (the scope).
 */
public interface SortedFactor<D extends Description> extends Factor<D> {

  /**
   * Return the hypothesis code corresponding to the affinity of given rank,
   * with rank 0 corresponding to the highest affinity in the factor.
   * @param rank position in sorted order
   * @return hypothesis code
   */
  int hypothesis(int rank);

  /**
   * Return the affinity in the factor arithmetic corresponding to the given rank,
   * with rank 0 corresponding to the highest affinity in the factor.
   * @param rank position in sorted order
   * @return affinity
   */
  double value(int rank);
}
