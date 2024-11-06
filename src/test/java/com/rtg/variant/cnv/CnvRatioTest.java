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
package com.rtg.variant.cnv;

import static com.rtg.util.StringUtils.LS;

import java.io.File;
import java.io.IOException;
import java.io.Reader;
import java.io.StringReader;
import java.util.HashMap;
import java.util.Map;

import com.rtg.launcher.OutputParams;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;
import com.rtg.variant.cnv.region.Region;
import com.rtg.variant.cnv.region.RegionUtils;

import junit.framework.TestCase;

/**
 */
public class CnvRatioTest extends TestCase {

  public void test() throws IOException {
    final File outputDir = FileUtils.createTempDir("cnvratio", "test");
    final Map<Integer, String> nameMap = new HashMap<>();
    nameMap.put(0, "simulatedSequence1");
    nameMap.put(1, "simulatedSequence2");
    nameMap.put(2, "simulatedSequence3");
    nameMap.put(3, "simulatedSequence4");
    nameMap.put(4, "simulatedSequence5");
    final String inStr = ""
      + "simulatedSequence4 2 3" + LS
      + "simulatedSequence4 4 5" + LS
      ;
    final Reader in = new StringReader(inStr);
    final Map<String, Region> regionMap = RegionUtils.regionsFromReader(in);

    final CnvRatio ratio = new CnvRatio(3.0, regionMap, nameMap, CnvProductParams.builder().outputParams(new OutputParams(outputDir, false)).bucketSize(5).create(), 10, false);
    try {
      ratio.exec(new int[][] {{3, 1, 2, 1, 0}, {3, 3, 2, 0, 0}, {3, 3, 1, 0, 0}, {3, 3, 2, 0, 0}, {9, 9, 9, 9, 9, 100, 100, 5, 5, 5, 5, 5, 5, 5, 5, 100, 100, 100}}, new int[][] {{1, 2, 1, 1, 0}, {3, 0, 2, 1, 0}, {0, 0, 4, 0, 0}, {3, 0, 2, 1, 0}, {9, 9, 9, 9, 9, 100, 100, 5, 5, 5, 5, 5, 5, 5, 5, 100, 100, 100}});
      final String bedFile = FileUtils.fileToString(new File(outputDir, "cnv.bed"));
      final String ratioFile = FileUtils.fileToString(new File(outputDir, "cnv.ratio"));
      final String bedBody = stripComments(bedFile);
      final String ratioBody = stripComments(ratioFile);
      final String bedComments = getComments(bedFile);
      final String ratioComments = getComments(ratioFile);
      assertEquals(StringUtils.convertLineEndings(FileHelper.resourceToString("com/rtg/variant/cnv/resources/cnv.cnv")), bedBody);
      assertEquals(StringUtils.convertLineEndings(FileHelper.resourceToString("com/rtg/variant/cnv/resources/cnv.ratio")), ratioBody);
      TestUtils.containsAll(bedComments, "#CL\t", "#Version ", ", cnv v",
          "#Seq\tstart\tend\tcnv\tmean\tstandard-deviation",
          "#Seq\tstart\tend\tnregion",
          "#Seq\tstart\tend\tgdelregion",
          "#Seq\tstart\tend\tmodel[0-5] mean");
      TestUtils.containsAll(ratioComments, "#CL\t", "#Version ", ", cnv v",
          "#Seq\tposition\tratio");
    } finally {
      ratio.closeOutputs();
      assertTrue(FileHelper.deleteAll(outputDir));
    }
  }

  public void testState() {
    TestUtils.testEnum(CnvRatio.State.class, "[IN, OUT]");
  }

  private String stripComments(String file) {
    final StringBuilder bd = new StringBuilder();
    final String[] lines = file.split(LS);
    for (String line : lines) {
      if (!line.startsWith("#")) {
        bd.append(line);
        bd.append(LS);
      }
    }
    return bd.toString();
  }

  private String getComments(String file) {
    final StringBuilder bd = new StringBuilder();
    final String[] lines = file.split(LS);
    for (String line : lines) {
      if (line.startsWith("#")) {
        bd.append(line);
        bd.append(LS);
      }
    }
    return bd.toString();
  }

}
