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
          final CombinedReport combined = new CombinedReport(Arrays.asList(dir), dir);
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
