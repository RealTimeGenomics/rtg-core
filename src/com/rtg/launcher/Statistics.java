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
package com.rtg.launcher;

import java.io.IOException;
import java.io.OutputStream;

/**
 */
public interface Statistics {

  /**
   * Prints the statistics in a human readable way to <code>reportStream</code>
   * @param reportStream stream to write to
   * @throws IOException if an exception occurs
   */
  void printStatistics(final OutputStream reportStream) throws IOException;

  /**
   * Generate an HTML report for this statistics type.
   * @throws IOException if an exception occurs during input or output
   */
  void generateReport() throws IOException;
}
