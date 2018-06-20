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
    try (TestDirectory temp = new TestDirectory()) {
      final File out = FileUtils.createTempDir("report", "out", temp);
      final File in = FileUtils.createTempDir("species", "in", temp);
      FileUtils.stringToFile(SIMILARITY_TSV, new File(in, "similarity.tsv"));

      final SimilarityReport report = new SimilarityReport();
      report.generateReport(ReportType.HTML, out, in);

      assertTrue(out.exists());
      final File index = new File(out, "index.html");
      assertTrue(index.exists());
      final File reportFilesDir = new File(out, "index_files");
      assertTrue(reportFilesDir.exists());
      assertTrue(new File(reportFilesDir, "rtg.css").exists());
      final File matrix = new File(reportFilesDir, "matrix.html");
      assertTrue(matrix.exists());

      final String matrixString = FileUtils.fileToString(matrix);
//      System.err.println(matrixString);
      TestUtils.containsAll(matrixString,
        "<tr><td>horse</td><td style=\"background-color: #00ff00\">",
        "<td style=\"background-color: #c03f00\">",
        "<td style=\"background-color: #40bf00\">",
        "<tr><td>cow</td><td style=\"background-color: #c03f00\">",
        "<td style=\"background-color: #00ff00\">",
        "<td style=\"background-color: #cc3300\">",
        "<tr><td>donkey</td><td style=\"background-color: #40bf00\">",
        "<td style=\"background-color: #cc3300\">",
        "<td style=\"background-color: #00ff00\">"
      );
    }
  }
}
