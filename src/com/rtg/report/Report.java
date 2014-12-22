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

package com.rtg.report;

import java.io.File;
import java.io.IOException;

/**
 */
public interface Report {
  /**
   * Creates a report for the specified type, reading from <code>inputDirs</code> and writing to <code>outputDir</code>.
   *
   * @param type type of report to generate
   * @param inputDirs directories to read data from
   * @param outputDir directory to write report to
   * @throws IOException if an I/O error occurs
   */
  void generateReport(ReportType type, File[] inputDirs, File outputDir) throws IOException;

  /**
   * Creates a report for the specified type, reading from <code>inputDir</code> and writing to <code>outputDir</code>.
   *
   * @param type type of report to generate
   * @param inputDir directory to read data from
   * @param outputDir directory to write report to
   * @throws IOException if an I/O error occurs
   */
  void generateReport(ReportType type, File inputDir, File outputDir) throws IOException;

}
