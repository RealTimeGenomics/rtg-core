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
import java.io.PrintStream;

import com.rtg.ngs.MapReportData;
import com.rtg.ngs.MapReportDataTest;
import com.rtg.sam.SamFilterParams;
import com.rtg.util.TestUtils;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;

/**
 */
public class MapReportTest extends AbstractReportTest {

  @Override
  Report getReport() {
    return new MapReport(SamFilterParams.builder().create());
  }

  public void testEnum() {
    TestUtils.testEnum(MapReport.PlotType.class, "[BOX, LINE, ORI, MORI, MAPC]");
  }
  public void testValidParams() throws Exception {
    try (TestDirectory testDir = new TestDirectory()) {
      final File out = FileHelper.createTempDirectory(testDir);
      final File in = FileHelper.createTempDirectory(testDir);
      final MapReportData.Merger data = new MapReportData.Merger();
      MapReportDataTest.createMapData(data);
      final MapReportData blended = data.blendReportData();
      try (final PrintStream inputFile = new PrintStream(FileUtils.createOutputStream(new File(in, MapReportData.MAP_REPORT_FILE_NAME)))) {
        blended.write(inputFile);
      }

      final MapReport report = new MapReport(SamFilterParams.builder().create());
      report.generateReport(ReportType.HTML, out, in);

      assertTrue(out.exists());

      final File index = new File(out, "index.html");
      assertTrue(index.exists());
      final String asString = FileUtils.fileToString(index);

      assertTrue(asString.contains("Mapped</td><td align=\"right\">10"));
      assertTrue(asString.contains("forward</td><td align=\"right\">10"));
      assertTrue(asString.contains("<strong>Command line: </strong>"));

      for (int i = 0; i < 4; ++i) {
        assertTrue(asString.contains(i + "</td><td align=\"right\">2"));
      }

      final File reportFilesDir = new File(out, "index_files");
      assertTrue(reportFilesDir.exists());
      assertTrue(new File(reportFilesDir, "rtg.css").exists());
    }
  }
}
