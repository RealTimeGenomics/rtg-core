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
import com.rtg.util.test.FileHelper;

/**
 */
public class ReadSimValidatorTest extends AbstractCliTest {

  @Override
  protected AbstractCli getCli() {
    return new ReadSimCli();
  }

  public void testCliValidator1() throws IOException, InvalidParamsException {
    Diagnostic.setLogStream();
    final File tempDir = FileUtils.createTempDir("readsimclitest", "checkcli");
    try {
      final File genomeDir = FileHelper.createTempDirectory();
      try {
        ReaderTestUtils.getReaderDNA(">seq1" + StringUtils.LS + "acgt", genomeDir, null).close();
        //final File xgenomeFile = new File(tempDir, "xgenome");
        final File reads = new File(tempDir, "reads");
        TestUtils.containsAll(checkHandleFlagsErr("-o", reads.getPath(), "-r", "36", "-n", "100"), "You must provide values for -t SDF --machine STRING");
        TestUtils.containsAll(checkHandleFlagsErr("-t", genomeDir.getPath(), "-r", "36", "-n", "100", "--machine", "illumina_se"), "You must provide a value for -o SDF");

        TestUtils.containsAll(checkHandleFlagsErr("-o", reads.getPath(), "-t", genomeDir.getPath(),
            "-R", "20", "-L", "20", "-r", "36", "-n", "100", "--machine", "illumina_pe"), "The flag --read-length is not permitted for this set of arguments");

        TestUtils.containsAll(checkHandleFlagsErr("-o", reads.getPath(), "-t", genomeDir.getPath(),
            "-R", "20", "-L", "20", "-n", "100", "--machine", "illumina_se"), "The flag --read-length is required");

        TestUtils.containsAll(checkHandleFlagsErr("-o", reads.getPath(), "-t", genomeDir.getPath(),
            "-n", "100", "--machine", "illumina_pe"), "The flag --left-read-length is required");
        TestUtils.containsAll(checkHandleFlagsErr("-o", reads.getPath(), "-t", genomeDir.getPath(),
            "--left-read-length", "30", "-n", "100", "--machine", "illumina_pe"), "The flag --right-read-length is required");

        TestUtils.containsAll(checkHandleFlagsErr("-t", genomeDir.getPath(), "-o", reads.getPath(), "-c", "0", "--machine", "illumina_se"), "Coverage should be positive");

        TestUtils.containsAll(checkHandleFlagsErr("-t", genomeDir.getPath(), "-o", reads.getPath(), "-c", "1000001", "--machine", "illumina_se"), "Coverage cannot be greater than 1000000.0");

        TestUtils.containsAll(checkHandleFlagsErr("-t", genomeDir.getPath(), "-o", reads.getPath(), "-n", "0", "--machine", "illumina_se"), "Number of reads should be greater than 0");

        TestUtils.containsAll(checkHandleFlagsErr("-t", genomeDir.getPath(), "-o", reads.getPath(), "-n", "1", "-r", "0", "--machine", "illumina_se"), "Read length is too small");

        TestUtils.containsAll(checkHandleFlagsErr("-t", genomeDir.getPath(), "-o", reads.getPath(), "-n", "1", "-r", "10", "--machine", "illumina_se", "--max-fragment-size", "10", "--min-fragment-size", "11"), "--max-fragment-size should not be smaller than --min-fragment-size");

        TestUtils.containsAll(checkHandleFlagsErr("-t", genomeDir.getPath(), "-o", reads.getPath(), "-n", "1", "-r", "11", "--machine", "illumina_se", "--max-fragment-size", "10", "--min-fragment-size", "10"), "Read length is too large for selected fragment size");

        TestUtils.containsAll(checkHandleFlagsErr("-t", genomeDir.getPath(), "-o", reads.getPath(), "-n", "1", "-r", "11", "--machine", "illumina_se", "-T", genomeDir.getPath()), "The --Xdiploid-input SDF cannot be the same as that given with --input");
      } finally {
        assertTrue(FileHelper.deleteAll(genomeDir));
      }

    } finally {
      assertTrue(FileHelper.deleteAll(tempDir));
    }
  }
}
