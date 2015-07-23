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

    final CnvRatio ratio = new CnvRatio(3.0, regionMap, nameMap, CnvProductParams.builder().outputParams(new OutputParams(outputDir, false, false)).bucketSize(5).create(), 10, false);
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
