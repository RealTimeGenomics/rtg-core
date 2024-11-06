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
public class ReadSimEvalParamsTest extends AbstractCliTest {

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
      ReadSimEvalParams.initFlags(flags);

      final String[] args = {"--output", rseDir.getPath(), "--reads", reads.getPath(),
          "--Xmismatch", "4", "-v", "3", "--Xgap-opening", "5", sam.getPath(),
          };

      assertTrue(flags.setFlags(args));

      final ReadSimEvalParams params = new ReadSimEvalParams(flags);

      assertEquals(3, params.variance());
      assertEquals(4, params.misMatchPenalty());
      assertEquals(5, params.gapOpeningPenalty());
      assertEquals(false, params.isPaired());
      assertEquals(false, params.dumpRecords());
      assertEquals(false, params.scoreHistograms());

      createPairedDir(reads);

      String[] args2 = {"--output", rseDir.getPath(), "-r", reads.getPath(), sam.getPath()};

      assertTrue(flags.setFlags(args2));

      final ReadSimEvalParams params2 = new ReadSimEvalParams(flags);

      assertEquals(0, params2.variance());
      assertEquals(1, params2.misMatchPenalty());
      assertEquals(2, params2.gapOpeningPenalty());
      assertEquals(true, params2.isPaired());
      assertEquals(false, params2.verbose());

      args2 = new String[] {"--output", rseDir.getPath(), "-r", reads.getPath(),
          sam.getPath(), "-m", "5", "--verbose"};

      assertTrue(flags.setFlags(args2));

      final ReadSimEvalParams params3 = new ReadSimEvalParams(flags);

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
    return new ReadSimEvalCli();
  }
}
