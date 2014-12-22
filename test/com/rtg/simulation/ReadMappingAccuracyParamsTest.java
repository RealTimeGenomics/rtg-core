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

package com.rtg.simulation;

import java.io.File;
import java.io.IOException;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.StringUtils;
import com.rtg.util.cli.CFlags;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;

/**
 */
public class ReadMappingAccuracyParamsTest extends AbstractCliTest {

  static final String READS = ">read0:1:1234:S1:I0:D0" + StringUtils.LS
     //0123456789012345678901234567890123456789012
    + "ACATGCTGCATGCATGCTGATGCTGCTGCTATAGTGATGAATC" + StringUtils.LS
    + ">read1:1:2345R:S0:I1:D0" + StringUtils.LS
    + "TGCTCATGTAGCTGAATTCTAGGACGCCATGCATGACTGACTG"  + StringUtils.LS
    + ">read2:1:56785:S4:I0:D0" + StringUtils.LS
    + "GTACTGCATCGATCGATCGTAGCTACGTAGCATGCATCGATGC" + StringUtils.LS;

  public void testAll() throws IOException, InvalidParamsException {
    Diagnostic.setLogStream();
    final File topLevel = FileUtils.createTempDir("readmappingaccuracy", "paramstest");
    try {
      final File reads = new File(topLevel, "reads");
      final File rseDir = new File(topLevel, "rse");
      ReaderTestUtils.getReaderDNA(READS, reads, null).close();
      final File sam = new File(topLevel, "samfile.sam");
      assertTrue(sam.createNewFile());

      final CFlags flags = new CFlags();
      ReadMappingAccuracyParams.initFlags(flags);

      final String[] args = {"--output", rseDir.getPath(), "--reads", reads.getPath(),
          "--Xmismatch", "4", "-v", "3", "--Xgap-opening", "5", sam.getPath(),
          };

      assertTrue(flags.setFlags(args));

      final ReadMappingAccuracyParams params = new ReadMappingAccuracyParams(flags);

      assertEquals(3, params.variance());
      assertEquals(4, params.misMatchPenalty());
      assertEquals(5, params.gapOpeningPenalty());
      assertEquals(false, params.isPaired());
      assertEquals(false, params.dumpRecords());
      assertEquals(false, params.scoreHistograms());

      createPairedDir(reads);

      String[] args2 = {"--output", rseDir.getPath(), "-r", reads.getPath(), sam.getPath()};

      assertTrue(flags.setFlags(args2));

      final ReadMappingAccuracyParams params2 = new ReadMappingAccuracyParams(flags);

      assertEquals(0, params2.variance());
      assertEquals(1, params2.misMatchPenalty());
      assertEquals(2, params2.gapOpeningPenalty());
      assertEquals(true, params2.isPaired());
      assertEquals(false, params2.verbose());

      args2 = new String[] {"--output", rseDir.getPath(), "-r", reads.getPath(),
          sam.getPath(), "-m", "5", "--verbose"};

      assertTrue(flags.setFlags(args2));

      final ReadMappingAccuracyParams params3 = new ReadMappingAccuracyParams(flags);

      assertEquals(5, params3.misMatchPenalty());
      assertTrue(params3.verbose());
    } finally {
      assertTrue(FileHelper.deleteAll(topLevel));
    }

  }

  private void createPairedDir(File toplevel) throws IOException {
    ReaderTestUtils.getReaderDNA(READS, new File(toplevel, "left"), null).close();
    ReaderTestUtils.getReaderDNA(READS, new File(toplevel, "right"), null).close();
  }

  public void testFlags() {
    checkHelp("SDF containing reads",
        "SAM/BAM format files",
        "variation allowed in start position",
        "exclude all mated",
        "exclude all unmated",
        "provide more detailed breakdown of stats");
  }

  @Override
  protected AbstractCli getCli() {
    return new ReadMappingAccuracy();
  }
}
