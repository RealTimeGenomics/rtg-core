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
import com.rtg.reader.ReaderTestUtils;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
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
    checkHelp("Simulate mutations in Complete Genomics reads.");
  }


  public void testCliValidator1() throws IOException, InvalidParamsException {
    Diagnostic.setLogStream();
    final File tempDir = FileUtils.createTempDir("readsimclitest", "checkcli");
    try {
      final File genomeDir = FileHelper.createTempDirectory();
      try {
        ReaderTestUtils.getReaderDNA(">seq1" + StringUtils.LS + "acgt", genomeDir, null).close();
        final File reads = new File(tempDir, "reads");
        TestUtils.containsAll(checkHandleFlagsErr("-o", reads.getPath(), "-n", "100"), "Usage: rtg cgsim [OPTION]... -t SDF -o SDF -c FLOAT", "You must provide a value for -t SDF");
        TestUtils.containsAll(checkHandleFlagsErr("-t", genomeDir.getPath(), "-n", "100"), "You must provide a value for -o SDF");
        TestUtils.containsAll(checkHandleFlagsErr("-t", genomeDir.getPath(), "-o", reads.getPath(), "-c", "0"), "Coverage should be positive");
        TestUtils.containsAll(checkHandleFlagsErr("-t", genomeDir.getPath(), "-o", reads.getPath(), "-n", "0"), "Number of reads should be greater than 0");
        assertEquals(MachineType.COMPLETE_GENOMICS, new CgSimCli().getMachineType());
        final CgSimCli.CgSimValidator validator = new CgSimCli.CgSimValidator();
        assertTrue(validator.checkMachines(null));
      } finally {
        assertTrue(FileHelper.deleteAll(genomeDir));
      }

    } finally {
      assertTrue(FileHelper.deleteAll(tempDir));
    }
  }
}
