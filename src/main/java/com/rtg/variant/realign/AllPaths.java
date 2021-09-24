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

package com.rtg.variant.realign;

import com.rtg.variant.util.arithmetic.PossibilityArithmetic;

/**
 */
public interface AllPaths {

  /**
   * Change the environment and recalculate all alignment paths.
   *
   * @param env the new read and template information.
   */
  void setEnv(final Environment env);

  /**
   * Get the overall probability of the read given the template.
   *
   * @return the natural log of the probability.
   */
  double totalScoreLn();

  /**
   * Get the overall probability of the read given the template.
   *
   * @return the probability as a normal double in the range from 0.0 to 1.0 inclusive.
   */
  double totalScore();

  /**
   * Get the overall probability of the read given the template.
   *
   * @return the probability as a possibility (see <code>arithmetic()</code>).
   */
  double total();

  /**
   * Check if the calculation has underflowed.
   * @return true iff the alignment paths calculation has underflowed.
   */
  boolean underflow();

  /**
   * @return the arithmetic object used for calculations.
   */
  PossibilityArithmetic arithmetic();
}
