/*
 * Copyright (c) 2016. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
 */

package com.rtg.variant.cnv.segment;

import java.io.File;
import java.io.IOException;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;
import com.rtg.util.TestUtils;

/**
 * Tests the corresponding class.
 */
public class SegmentCliTest extends AbstractCliTest {

  @Override
  protected AbstractCli getCli() {
    return new SegmentCli();
  }

  public void testHelp() {
    checkHelp("Segments depth of coverage data to identify regions of consistent copy number.",
      "BED file supplying per-region coverage data for the sample",
      "BED file supplying per-region coverage data for control sample",
      "directory for output",
      "SDF containing reference genome",
      "weighting factor for inter-segment distances during energy scoring",
      "weighting factor for intra-segment distances during energy scoring",
      "segmentation sensitivity"
      );
  }

  public void testErrors() throws IOException {
    final File emptyFile = File.createTempFile("test", ".vcf");
    try {
      String res = checkHandleFlagsErr();
      final String exp = getCFlags().getUsageHeader();
      assertTrue(res.contains(exp));
      assertTrue(res.contains("Error: You must provide values for --case FILE -o DIR -t SDF"));
      res = checkHandleFlagsErr("-o", "test-foo-out", "-t", "test-sdf", "--case", emptyFile.getPath());
      assertTrue(res.contains(exp));
      TestUtils.containsAll(res, "Error: One of --Xcolumn or --control or --panel must be set");
      res = checkHandleFlagsErr("-o", "test-foo-out", "-t", "test-sdf", "--case", emptyFile.getPath(), "--control", emptyFile.getPath(), "--Xlimit", "0");
      assertTrue(res.contains(exp));
      assertTrue(res.contains("Error: The value for --Xlimit must be at least 1"));
    } finally {
      assertTrue(emptyFile.delete());
    }
  }

}
