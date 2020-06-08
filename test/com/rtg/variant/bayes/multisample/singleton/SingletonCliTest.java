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
