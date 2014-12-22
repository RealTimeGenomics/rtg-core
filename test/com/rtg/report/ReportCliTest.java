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

import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.IOException;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;
import com.rtg.util.StringUtils;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;

/**
 */
public class ReportCliTest extends AbstractCliTest {

  @Override
  public void setUp() throws IOException {
    super.setUp();
  }

  @Override
  public void tearDown() throws IOException {
    super.tearDown();
  }

  @Override
  protected AbstractCli getCli() {
    return new ReportCli();
  }

  public void testInValidFlags() {
    checkHandleFlagsErr();
  }

  public void testValidParams() throws Exception {
    final File temp = FileUtils.createTempDir("test", "temp");
    try {
      final File out = new File(temp, "out");
      final File in = FileHelper.createTempDirectory(temp);
      final String content = "#fraction-individuals\tfraction-DNA\tdepth\tbreadth\tname" + StringUtils.LS
          + "0.312221\t0.276982\t0.552580\t0.001676\tgi|15805042|ref|NC_001263.1|" + StringUtils.LS
          + "0.126791\t0.096501\t0.537990\t0.002665\tgi|77358697|ref|NC_003112.2|" + StringUtils.LS;
      final ByteArrayInputStream bais = new ByteArrayInputStream(content.getBytes());
      FileHelper.streamToFile(bais, new File(in, "species.tsv"));

      checkMainInitOk("-o", out.getPath(), "-t", "html", "-m", "species", in.getPath());
      assertTrue(out.exists());
      assertTrue(new File(out, "done").exists());
      assertTrue(new File(out, "index.html").exists());

      final File reportFilesDir = new File(out, "index_files");
      assertTrue(reportFilesDir.exists());
      assertTrue(new File(reportFilesDir, "rtg.css").exists());
      assertTrue(new File(reportFilesDir, "rtg_logo.png").exists());
    } finally {
      assertTrue(FileHelper.deleteAll(temp));
    }
  }

  public void testInitParams() {
    checkHelp("report [OPTION]... -m STRING -o DIR DIR+"
        , "Creates a report for outputs from RTG modules."
        , "-o", "--output=DIR", "directory to write report to"
        , "DIR+", "RTG module output directory. Must be specified 1 or more times"
        , "-m", "--module=STRING", "module to generate report for (Must be one of [species, similarity, map])"
        , "-t", "--report-type=STRING", "type of report to generate (Must be one of [html, pdf]) (Default is html)"
        );
  }
}
