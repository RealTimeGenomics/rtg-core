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
package com.rtg.simulation.reads;

import java.io.File;
import java.io.IOException;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;
import com.rtg.launcher.MainResult;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.machine.MachineType;
import com.rtg.util.test.FileHelper;

/**
 * Test class
 */
public class CgSimCliTest extends AbstractCliTest {

  @Override
  protected AbstractCli getCli() {
    return new CgSimCli();
  }

  public void testFlags() {
    checkHelp("Simulate Complete Genomics Inc sequencing reads.");
  }


  public void testCliValidator1() throws IOException, InvalidParamsException {

    try (final TestDirectory tempDir = new TestDirectory("cgsimclitest")) {
      final File genomeDir = FileHelper.createTempDirectory();
      try {
        ReaderTestUtils.getReaderDNA(">seq1" + StringUtils.LS + "acgt", genomeDir, null).close();
        final File reads = new File(tempDir, "reads");
        TestUtils.containsAllUnwrapped(checkHandleFlagsErr("-o", reads.getPath(), "-n", "100"), "Usage: rtg cgsim [OPTION]... -V INT -t SDF -o SDF -c FLOAT", "You must provide values for -V INT -t SDF");
        TestUtils.containsAllUnwrapped(checkHandleFlagsErr("-t", genomeDir.getPath(), "-n", "100"), "You must provide values for -V INT -o SDF");
        TestUtils.containsAllUnwrapped(checkHandleFlagsErr("-t", genomeDir.getPath(), "-o", reads.getPath(), "-c", "0", "-V", "1"), "Coverage should be positive");
        TestUtils.containsAllUnwrapped(checkHandleFlagsErr("-t", genomeDir.getPath(), "-o", reads.getPath(), "-n", "0", "-V", "1"), "Number of reads should be greater than 0");
        TestUtils.containsAllUnwrapped(checkHandleFlagsErr("-t", genomeDir.getPath(), "-o", reads.getPath(), "-n", "1", "--cg-read-version", "0"), "Version must be 1 or 2");
        TestUtils.containsAllUnwrapped(checkHandleFlagsErr("-t", genomeDir.getPath(), "-o", reads.getPath(), "-n", "1", "--cg-read-version", "3"), "Version must be 1 or 2");

        final CgSimCli cli = (CgSimCli) getCli();
        MainResult.run(cli, "-t", genomeDir.getPath(), "-o", reads.getPath(), "-n", "1", "-V", "1");
        assertEquals(MachineType.COMPLETE_GENOMICS, cli.getMachineType());
      } finally {
        assertTrue(FileHelper.deleteAll(genomeDir));
      }
    }
  }
}
