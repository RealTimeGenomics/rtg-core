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
 * Store the measure for a sample's probability distribution
 */
public interface GenotypeMeasure {
  /**
   * @return the arithmetic  used by this measure
   */
  PossibilityArithmetic arithmetic();

  /**
   *
   * @param hypothesis the hypothesis of interest
   * @return the un-normalised measure for that hypothesis
   */
  double measure(int hypothesis);

  /**
   * @return the number of hypotheses
   */
  int size();

  /**
   *
   * @return the posterior score of the best hypothesis
   */
  double bestPosterior();

  /**
   * @return the posterior of the call not being the reference hypothesis
   */
  double nonIdentityPosterior();


  /**
   *
   * @return the id of the best hypothesis
   */
  int best();

  /**
   * @return the id of the best hypothesis
   */
  int reference();

  /**
   * @return the hypothesis space this measure is over
   */
  Hypotheses<?> hypotheses();
}
