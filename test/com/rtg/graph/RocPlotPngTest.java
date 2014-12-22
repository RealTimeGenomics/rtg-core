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
package com.rtg.graph;

import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.Collections;

import javax.imageio.ImageIO;

import com.rtg.util.StringUtils;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.TestDirectory;

import junit.framework.TestCase;

/**
 */
public class RocPlotPngTest extends TestCase {

  private static final String ROC = ""
          + "#total baseline variants: 3092754" + StringUtils.LS
          + "#score\ttrue_positives\tfalse_positives" + StringUtils.LS
          + "3.300\t0.000\t15" + StringUtils.LS
          + "2.261\t70000.000\t137" + StringUtils.LS
          + "1.226\t180000.000\t516" + StringUtils.LS
          + "0.700\t406000.000\t11337" + StringUtils.LS
          + "0.533\t1971000.000\t1446920" + StringUtils.LS
          + "0.333\t2071000.000\t1646920" + StringUtils.LS
          + "0.200\t2995295.000\t1864591" + StringUtils.LS;

  public void test() throws IOException {
    try (final TestDirectory dir = new TestDirectory()) {
      final File roc = FileUtils.stringToFile(ROC, new File(dir, "roc.tsv"));
      final File png = new File(dir, "PNG.png");
      RocPlotPng.rocPngImage(Collections.singletonList(roc), Collections.singletonList("LINE"), "a title", true, 3, png);
      final BufferedImage buf = ImageIO.read(png);
      assertEquals(800, buf.getWidth());
      assertEquals(600, buf.getHeight());
    }
  }
}
