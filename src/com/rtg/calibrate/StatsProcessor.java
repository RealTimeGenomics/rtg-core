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

package com.rtg.calibrate;

/**
 * Interface for call-back functions for <code>Calibrator.processStats</code>.
 *
 */
public interface StatsProcessor {

  /**
   * Call-back with covariate values with indexes corresponding to covariate index.
   * @param covariateValues the current indices into the hypercube.
   * @param stats the total counts at that position.  Null means all the counts are zero.
   */
  void process(int[] covariateValues, CalibrationStats stats);
}
