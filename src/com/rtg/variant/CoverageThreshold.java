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
package com.rtg.variant;

/**
 */
public interface CoverageThreshold {

  /**
   * Get the single-sample coverage threshold for the requested sequence
   * @param sequenceName name of the sequence to threshold
   * @return the maximum coverage allowed for the sequence for a single sample
   */
  int thresholdSingle(String sequenceName);

  /**
   * Get the all-sample coverage threshold for the requested sequence
   * @param sequenceName name of the sequence to threshold
   * @return the maximum coverage allowed for the sequence across all samples
   */
  int thresholdTotal(String sequenceName);

}
