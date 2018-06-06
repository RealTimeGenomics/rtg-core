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

import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;

/**
 */
public class SimilarityReportTest extends AbstractReportTest {
  @Override
  Report getReport() {
    return new SimilarityReport();
  }
  static final String SIMILARITY_TSV = (" horse cow donkey" + StringUtils.LS
      + "horse 20 5 15" + StringUtils.LS
      + "cow 5 20 4" + StringUtils.LS
      + "donkey 15 4 20" + StringUtils.LS).replaceAll(" ", "\t");

  public void testValidParams() throws Exception {
    try (TestDirectory testDir = new TestDirectory()) {
      final File out = FileHelper.createTempDirectory(testDir);
      final File in = FileHelper.createTempDirectory(testDir);
      FileUtils.stringToFile(SIMILARITY_TSV, new File(in, "similarity.tsv"));

      final SimilarityReport report = new SimilarityReport();
      report.generateReport(ReportType.HTML, out, in);

      assertTrue(out.exists());
//      System.err.println(Arrays.toString(out.listFiles()));
      final File index = new File(out, "index.html");
      assertTrue(index.exists());
      final File reportFilesDir = new File(out, "index_files");
      assertTrue(reportFilesDir.exists());
      assertTrue(new File(reportFilesDir, "rtg.css").exists());
      assertTrue(new File(reportFilesDir, "rtg_logo.png").exists());
      final File matrix = new File(reportFilesDir, "matrix.html");
      assertTrue(matrix.exists());

      final String matrixString = FileUtils.fileToString(matrix);
//      System.err.println(matrixString);
      TestUtils.containsAll(matrixString
          , "<tr><td>horse</td><td style=\"background-color: #00ff00\"><img src=\"blank.png\" width=\"100\" height=\"15\" title=\"20\" /></td><td style=\"background-color: #c03f00\"><img src=\"blank.png\" width=\"100\" height=\"15\" title=\"5\" /></td><td style=\"background-color: #40bf00\"><img src=\"blank.png\" width=\"100\" height=\"15\" title=\"15\" /></td></tr>"
          , "<tr><td>cow</td><td style=\"background-color: #c03f00\"><img src=\"blank.png\" width=\"100\" height=\"15\" title=\"5\" /></td><td style=\"background-color: #00ff00\"><img src=\"blank.png\" width=\"100\" height=\"15\" title=\"20\" /></td><td style=\"background-color: #cc3300\"><img src=\"blank.png\" width=\"100\" height=\"15\" title=\"4\" /></td></tr>"
          , "<tr><td>donkey</td><td style=\"background-color: #40bf00\"><img src=\"blank.png\" width=\"100\" height=\"15\" title=\"15\" /></td><td style=\"background-color: #cc3300\"><img src=\"blank.png\" width=\"100\" height=\"15\" title=\"4\" /></td><td style=\"background-color: #00ff00\"><img src=\"blank.png\" width=\"100\" height=\"15\" title=\"20\" /></td></tr>"
      );
    }
  }
}
