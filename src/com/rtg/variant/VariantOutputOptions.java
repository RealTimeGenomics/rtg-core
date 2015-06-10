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
public interface VariantOutputOptions {

  /**
   * Return the level of calls that should be output.
   * @return calling level
   */
  VariantOutputLevel callLevel();

  /**
   * Whether to attempt calls in complex regions.
   * If true then no calls are attempted in complex regions but a single line with type "x" is output.
   * @return a boolean.
   */
  boolean noComplexCalls();

  /**
   * Only output if coverage is below this threshold.
   * @return coverage threshold
   */
  CoverageThreshold maxCoverageFilter();

  /**
   * Get the threshold for <code>NH/IH</code> ratio filtering.
   * @return the threshold for <code>NH/IH</code> ratio filtering (null if no filtering to be done).
   */
  Double maxAmbiguity();

}
