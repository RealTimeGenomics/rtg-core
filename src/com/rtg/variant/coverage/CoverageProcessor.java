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

package com.rtg.variant.coverage;

import java.io.IOException;

/**
 *         Date: 16/02/12
 *         Time: 1:45 PM
 */
public interface CoverageProcessor {
  /**
   * Performs any required initialisation
   * @throws IOException because it might do some output
   */
  void init() throws IOException;

  /**
   * @param name sequence name
   * @param position position on sequence
   * @param ih1 count of <code>IH=1</code> records
   * @param ihgt1 count of <code>IH>1</code> records
   * @param coverage coverage at position
   * @throws IOException if an IO error occurs
   */
  void finalCoveragePosition(String name, int position, int ih1, int ihgt1, double coverage) throws IOException;

  /**
   *
   * @param name sequence name
   * @param start start position
   * @param end end position
   * @param coverage coverage for the region
   * @throws IOException because it might do some output
   */
  void finalCoverageRegion(String name, int start, int end, int coverage) throws IOException;

  /**
   * Clean up any open resources
   * @throws IOException because it might do some output
   */
  void close() throws IOException;
}
