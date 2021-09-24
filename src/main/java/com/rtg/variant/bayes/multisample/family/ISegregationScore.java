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

package com.rtg.variant.bayes.multisample.family;

/**
 */
public interface ISegregationScore {

  /**
   * Increment the count of a particular hypothesis for a child
   * @param child the child hypothesis code
   * @param diploid true if the child is diploid, false otherwise
   */
  void increment(int child, boolean diploid);

  /**
   * Calculate the segregation log probability score
   * @return the segregation log probability score
   */
  double lnProbability();

}
