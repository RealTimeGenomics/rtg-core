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
