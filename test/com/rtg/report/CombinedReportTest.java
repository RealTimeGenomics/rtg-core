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
import java.util.Arrays;

import com.rtg.util.io.FileUtils;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;

/**
 */
public class CombinedReportTest extends TestCase {
  public void testCommandLine() throws IOException {
    try (TestDirectory dir = new TestDirectory()) {
      final File file = new File(dir, "mapx.log");
      FileHelper.resourceToFile("com/rtg/report/resources/mapx.log.txt", file);
      assertEquals("rtg mapx --output /tmp/metagenomics/mapx1 --template /tmp/protein_actual --max-alignment-score 10% --max-top-results 10 --input /tmp/metagenomics/mapf/unmapped.sdf/left", CombinedReport.getCommandLine(file));
    }
  }

  public void testLogDetect() throws IOException {
    final String[] types = {"mapf", "mapx", "species", "map"};
    for (String type : types) {
      if (!type.equals("species")) {
        try (TestDirectory dir = new TestDirectory()) {
          final File log = new File(dir, type + ".log");
          if (!log.createNewFile()) {
            throw new IOException();
          }
          final CombinedReport combined = new CombinedReport(Arrays.<File>asList(dir), dir);
//          System.err.println(type);
          combined.makeReport();
          final String s = FileUtils.fileToString(new File(dir, "index.html"));
          for (String checkType : types) {
            if (checkType.equals(type)) {
              final String expected = "class=\"checked\"> RTG " + checkType;
              assertTrue("<" + expected + "> was not contained in <" + s + ">", s.contains(expected));
            } else {
              final String expected = "class=\"crossed\"> RTG " + checkType;
              assertTrue("<" + expected + "> was not contained in <" + s + ">", s.contains(expected));
            }
          }
        }
      }
    }
  }

}
