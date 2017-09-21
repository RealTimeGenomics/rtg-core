/*
 * Copyright (c) 2017. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
 */

package com.rtg.variant.bayes.multisample.cancer;

/**
 * Provides a way to get the priors for cancer calling.
 */
public interface SomaticPriorsFactory {

  /**
   * Return the somatic Q prior matrix for the specified mutation rate.
   * @param mu probability of a somatic mutation.
   * @return probabilities of somatic transitions between possibly diploid hypotheses.
   */
  double[][] somaticQ(final double mu);

}
