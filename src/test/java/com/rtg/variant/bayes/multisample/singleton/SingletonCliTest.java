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

package com.rtg.variant.bayes.multisample.singleton;

import static com.rtg.util.StringUtils.LS;

import java.io.File;
import java.io.IOException;

import com.rtg.launcher.AbstractCli;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.tabix.TabixIndexer;
import com.rtg.tabix.UnindexableDataException;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.TestUtils;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;
import com.rtg.variant.bayes.multisample.AbstractCallerCliTest;
import com.rtg.variant.bayes.multisample.AbstractMultisampleCli;

/**
 */
public class SingletonCliTest extends AbstractCallerCliTest {

  @Override
  protected AbstractCli getCli() {
    return new SingletonCli();
  }

  @Override
  protected String getModuleName() {
    return "snp";
  }

  private static final String REF = ""
      + ">simulatedSequence1" + LS
      + "CATCATAACTGTCGACCAACACGGAGTCCACATCCCTTATCGGGACTCATCGGGTGGGACACTTGAGTCCGACCTGCGGATTAACGTATACAGTCGGCTG" + LS
      + ">simulatedSequence2" + LS
      + "AAAACGAGTGTATGAGGAAATGCGACAGCTACCCCCACCCGATTTAGCTGGCGGTTGCCGCCCTACGAGAAGATTTCTGCGCACAACCTTCGTCTCATTG" + LS;

  public void testCoverageLoading() throws InvalidParamsException, IOException, UnindexableDataException {
    try (final TestDirectory dir = new TestDirectory()) {
      final File ref = new File(dir, "ref");
      ReaderTestUtils.getReaderDNA(REF, ref, null).close();
      final File mated = new File(dir, "mated.sam.gz");
      final File unmated = new File(dir, "unmated.sam.gz");
      FileHelper.resourceToFile("com/rtg/variant/resources/coverage_mated.sam.gz", mated);
      FileHelper.resourceToFile("com/rtg/variant/resources/coverage_mated.sam.gz.calibration", new File(dir, "mated.sam.gz.calibration"));
      new TabixIndexer(mated, new File(dir, "mated.sam.gz.tbi")).saveSamIndex();
      FileHelper.resourceToFile("com/rtg/variant/resources/coverage_unmated.sam.gz", unmated);
      FileHelper.resourceToFile("com/rtg/variant/resources/coverage_unmated.sam.gz.calibration", new File(dir, "unmated.sam.gz.calibration"));
      new TabixIndexer(unmated, new File(dir, "unmated.sam.gz.tbi")).saveSamIndex();

      final File out = new File(dir, "calls");
      checkMainInitWarn("-t", ref.getPath(), "-o", out.getPath(), mated.getPath(), unmated.getPath(),
        "--filter-depth-multiplier", "3.0");
      final String log = FileUtils.fileToString(new File(out, "snp.log"));
      TestUtils.containsAll(log, "max_coverage_filter=avgSeqCov*3.0", //"max_coverage_threshold=avgSeqCov+(sqrt(avgSeqCov)*3.0)",
        "Sequence simulatedSequence2 filter on maximum per-sample coverage is 7",
        "Sequence simulatedSequence1 filter on maximum per-sample coverage is 10" //9"
      );
    }
  }

  public void testBug1508() throws InvalidParamsException, IOException, UnindexableDataException {
    // tests an issue in the circular buffer
    try (final TestDirectory dir = new TestDirectory()) {
      final String refFasta = FileHelper.resourceToString("com/rtg/variant/resources/bug1508_ref.fasta");
      final File ref = new File(dir, "ref");
      ReaderTestUtils.getReaderDNA(refFasta, ref, null).close();
      final File reads = new File(dir, "reads.sam.gz");
      FileHelper.resourceToFile("com/rtg/variant/resources/bug1508_reads.sam.gz", reads);
      new TabixIndexer(reads, new File(dir, "reads.sam.gz.tbi")).saveSamIndex();

      final File out = new File(dir, "calls");
      checkMainInitOk("-t", ref.getPath(), "-o", out.getPath(), reads.getPath(), "--" + AbstractMultisampleCli.NO_CALIBRATION);
      final String log = FileUtils.fileToString(new File(out, "snp.log"));
      TestUtils.containsAll(log, "16 alignments processed",
          "Finished successfully"
          );
    }
  }

  public void testBug1675NoRefCallWithAllMode() throws Exception {
    try (final TestDirectory dir = new TestDirectory()) {
      final String refFasta = FileHelper.resourceToString("com/rtg/variant/resources/bug1675_ref.fasta");
      final File ref = new File(dir, "ref");
      ReaderTestUtils.getReaderDNA(refFasta, ref, null).close();
      final File reads = new File(dir, "reads.sam.gz");
      FileHelper.resourceToFile("com/rtg/variant/resources/bug1675_reads.sam.gz", reads);
      new TabixIndexer(reads, new File(dir, "reads.sam.gz.tbi")).saveSamIndex();
      final File out = new File(dir, "calls");
      checkMainInitOk("-t", ref.getPath(), "-o", out.getPath(), reads.getPath(), "--" + AbstractMultisampleCli.NO_CALIBRATION, "--all", "--max-coverage", "25", "--region", "4:991+1");
      final String log = FileUtils.fileToString(new File(out, "snp.log"));
      TestUtils.containsAll(log, "160 alignments processed", "Finished successfully");
    }
  }
}
