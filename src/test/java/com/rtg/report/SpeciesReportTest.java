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

import java.io.ByteArrayInputStream;
import java.io.File;

import com.rtg.util.StringUtils;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;

/**
 */
public class SpeciesReportTest extends AbstractReportTest {

  @Override
  Report getReport() {
    return new SpeciesReport();
  }

  public void testValidParams() throws Exception {
    try (final TestDirectory temp = new TestDirectory("speciesreport")) {
      final File out = FileUtils.createTempDir("report", "out", temp);
      final File in = FileUtils.createTempDir("species", "in", temp);
      final String content = "#fraction-individuals\tfraction-DNA\tdepth\tbreadth\tname" + StringUtils.LS
          + "0.312221\t0.276982\t0.552580\t0.001676\tgi|15805042|ref|NC_001263.1|" + StringUtils.LS
          + "0.126791\t0.096501\t0.537990\t0.002665\tgi|77358697|ref|NC_003112.2|" + StringUtils.LS
          + "# some other junk" + StringUtils.LS;
      final ByteArrayInputStream bais = new ByteArrayInputStream(content.getBytes());
      FileHelper.streamToFile(bais, new File(in, "species.tsv"));

      final SpeciesReport report = new SpeciesReport();
      report.generateReport(ReportType.HTML, out, in);

      assertTrue(out.exists());
      assertTrue(new File(out, "index.html").exists());

      final File reportFilesDir = new File(out, "index_files");
      assertTrue(reportFilesDir.exists());
      assertTrue(new File(reportFilesDir, "rtg.css").exists());
      assertTrue(new File(reportFilesDir, "table.css").exists());
      assertTrue(new File(reportFilesDir, "table.js").exists());
    }
  }

  public void testTwoPointOh() throws Exception {
    try (final TestDirectory temp = new TestDirectory("speciesreport")) {
      final File out = FileUtils.createTempDir("report", "out", temp);
      final File in = FileUtils.createTempDir("species", "in", temp);
      final String content =  "#Version vPOST-2.7 build <not found> (<not found>), species v2.0\n"
                                    + "#abundance\tabundance-low\tabundance-high\tDNA-fraction\tDNA-fraction-low\tDNA-fraction-high\tconfidence\tcoverage-depth\tcoverage-breadth\tmapped-reads\ttaxon-id\tparent-id\trank\ttaxonomy-name\n"
                                    + "1.000\t0.9989\t1.000\t0.8019\t0.8010\t0.8028\t4154\t0.06472\t0.05506\t1856269.00\t131567\t1\tno rank\tcellular organisms\n"
                                    + "0.8493\t0.8483\t0.8503\t0.7268\t0.7260\t0.7277\t3949\t0.07165\t0.06286\t1680562.14\t2\t131567\tsuperkingdom\tBacteria\n"
                                    + "0.8408\t0.8398\t0.8417\t0.7175\t0.7167\t0.7184\t3932\t0.1411\t0.1164\t1658790.43\t1239\t2\tphylum\tFirmicutes\n"
                                    + "0.8234\t0.8224\t0.8244\t0.6895\t0.6887\t0.6903\t3852\t0.1429\t0.1201\t1593218.65\t91061\t1239\tclass\tBacilli" + StringUtils.LS;

      final ByteArrayInputStream bais = new ByteArrayInputStream(content.getBytes());
      FileHelper.streamToFile(bais, new File(in, "species.tsv"));

      final SpeciesReport report = new SpeciesReport();
      report.generateReport(ReportType.HTML, out, in);

      assertTrue(out.exists());
      assertTrue(new File(out, "index.html").exists());

      final File reportFilesDir = new File(out, "index_files");
      assertTrue(reportFilesDir.exists());
      assertTrue(new File(reportFilesDir, "rtg.css").exists());
      assertTrue(new File(reportFilesDir, "table.css").exists());
      assertTrue(new File(reportFilesDir, "table.js").exists());
    }
  }
}
