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
import com.rtg.launcher.MainResult;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.util.StringUtils;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;

/**
 * Tests the corresponding class.
 */
public class CnvPonBuildCliTest extends AbstractCliTest {

  @Override
  protected AbstractCli getCli() {
    return new CnvPonBuildCli();
  }

  public void testHelp() {
    checkHelp("Construct a normalized coverage sample from a panel of coverage outputs.",
      "BED output file",
      "SDF containing reference genome",
      "coverage BED file. Must be specified 1 or more times",
      "do not gzip the output",
      "do not produce indexes for output files",
      "Utility",
      "File Input/Output"
      );
  }

  public void testErrors() throws IOException {
    final String res = checkHandleFlagsErr();
    final String exp = getCFlags().getUsageHeader();
    assertTrue(res.contains(exp));
    assertTrue(res.contains("Error: You must provide values for -o FILE -t SDF FILE+"));
  }

  public void testValid() throws IOException {
    try (TestDirectory dir = new TestDirectory()) {
      final File coverage = FileHelper.resourceToFile("com/rtg/variant/cnv/segment/resources/coverage.bed", new File(dir, "coverage.bed"));
      final File output = new File(dir, "panel.bed");
      // Slightly G rich reference so G+C normalization actually does something
      final String ref = ">1" + System.getProperty("line.separator") + StringUtils.repeat("ACGGT", 25000 / "ACGGT".length());
      final File sdf = new File(dir, "sdf");
      ReaderTestUtils.getDNADir(ref, sdf);
      final MainResult result = checkMainInit("-Z", "-o", output.getPath(), "-t", sdf.getPath(), coverage.getPath());
      final String actual = StringUtils.grepMinusV(FileHelper.fileToString(output), "^#");
      mNano.check("expected.panel.bed", actual);
      mNano.check("expected.out.txt", result.out().replaceAll("/[^ ]*/", ""));
    }
  }

}
