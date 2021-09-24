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

package com.rtg.metagenomics;

/**
 * Detect termination in solve loops.
 * Needs to be done different ways if normal convergence or P-value convergence.
 */
interface Terminator {

  /**
   * @param estimate of the total accumulated L-value from here to infinity.
   * @param current L value.
   * @return true iff should stop solve iterations.
   */
  boolean terminate(double estimate, double current);
}
