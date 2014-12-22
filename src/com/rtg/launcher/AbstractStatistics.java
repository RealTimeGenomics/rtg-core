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

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;

import com.rtg.util.HtmlReportHelper;
import com.rtg.util.StringUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;

/**
 * Abstract statistics class with common code for outputting to log and stream
 */
public abstract class AbstractStatistics implements Statistics {

  static {
    System.setProperty("java.awt.headless", "true");
  }

  private final File mOutputDirectory;

  /**
   * @param outputDirectory The base output directory to generate statistics and reports in. May be null if no statistics or reports are to be generated.
   */
  public AbstractStatistics(File outputDirectory) {
    mOutputDirectory = outputDirectory;
  }

  protected HtmlReportHelper getReportHelper() {
    if (mOutputDirectory != null) {
      return new HtmlReportHelper(mOutputDirectory, "index");
    }
    return null;
  }

  /**
   * Method to get the statistics to print or log.
   * May return null if nothing is to be printed or logged.
   * @return the statistics string or null.
   */
  protected abstract String getStatistics();

  @Override
  public void printStatistics(OutputStream reportStream) throws IOException {
    final String stats = getStatistics();
        if (stats != null) {
          if (reportStream != null) {
            //print to output stream (std out often)
            reportStream.write(stats.getBytes());

            //now log to diagnostic
            final String[] lines = stats.split(StringUtils.LS);
            for (final String line : lines) {
              Diagnostic.userLog(line);
            }
      }
      //write to summary file if we have an output directory
      if (mOutputDirectory != null && mOutputDirectory.isDirectory()) {
        try (OutputStream summaryFile = FileUtils.createOutputStream(new File(mOutputDirectory, CommonFlags.SUMMARY_FILE), false)) {
          summaryFile.write(stats.getBytes());
        }
      }
    }
  }
}
