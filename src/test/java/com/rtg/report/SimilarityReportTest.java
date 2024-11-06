/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
