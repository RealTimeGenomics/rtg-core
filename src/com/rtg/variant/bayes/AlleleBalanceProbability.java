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
 * Possible implementations of the allele balance probability factor
 */
public interface AlleleBalanceProbability {

  /**
   * @param i current hypothesis
   * @param hypotheses list of hypotheses
   * @param statistics statistics at current location
   * @return the log scaled probability of the observed allele balance given the current hypothesis
   */
  double alleleBalanceLn(int i, Hypotheses<?> hypotheses, Statistics<?> statistics);
}
