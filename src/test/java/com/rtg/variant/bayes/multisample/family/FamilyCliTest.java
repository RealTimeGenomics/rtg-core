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
package com.rtg.variant.bayes.multisample.family;

import static com.rtg.sam.SharedSamConstants.SAM_BODY;
import static com.rtg.sam.SharedSamConstants.SAM_FAMILY;
import static com.rtg.util.StringUtils.LS;

import java.io.File;

import com.rtg.calibrate.RecalibrateCli;
import com.rtg.launcher.AbstractParamsCliTest;
import com.rtg.launcher.MainResult;
import com.rtg.launcher.ParamsCli;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.reference.Ploidy;
import com.rtg.reference.ReferenceTextBuilder;
import com.rtg.reference.Sex;
import com.rtg.tabix.TabixIndexer;
import com.rtg.util.TestUtils;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;
import com.rtg.variant.VariantParams;
import com.rtg.variant.bayes.multisample.AbstractMultisampleCli;

/**
 */
public class FamilyCliTest extends AbstractParamsCliTest<VariantParams> {

  @Override
  protected ParamsCli<VariantParams> getParamsCli() {
    return new FamilyCli();
  }

  private static final String EXP_F1 = "Error: You must provide values for -o DIR -t SDF" + LS;

  public void testErrorF1() {
    assertTrue(checkHandleFlagsErr().contains(EXP_F1));
  }

  public void testInitParams() {
    checkHelp("family [OPTION]... -o DIR -t SDF --father STRING --mother STRING FILE+",
        "sample identifier used in read groups for a daughter",
        "threshold for ambiguity"
    );
  }

  /**
   * Sam file header
   */
  public static final String SAMHEADER = ""
    + "@HD\tVN:1.0\tSO:coordinate" + LS
    + "@SQ\tSN:g1\tLN:20" + LS
    + "@RG\tID:RG1\tSM:TEST\tPL:ILLUMINA" + LS
    ;

  private static final String FAMILY_PED = ""
    + "0\tfather\t0\t0\t1\t0" + LS
    + "0\tmother\t0\t0\t2\t0" + LS
    + "0\tchild\tfather\tmother\t1\t0" + LS
    ;

  public void testValidator() throws Exception {
    final File tmpDir = FileHelper.createTempDirectory();
    final File tmpFile = FileUtils.stringToFile("original-derived TEST cancer contamination=0.13", FileHelper.createTempFile());
    final File tmpFile2 = FileUtils.stringToFile(FAMILY_PED, FileHelper.createTempFile());
    try {
      final File in = new File(tmpDir, "alignments.sam.gz");
      FileHelper.stringToGzFile(SAM_FAMILY, in);
      new TabixIndexer(in, new File(tmpDir, "alignments.sam.gz.tbi")).saveSamIndex();

      final File template = new File(tmpDir, "template");
      ReaderTestUtils.getDNADir(">g1\nacgtacgtacgtacgtacgt", template);
      final File outDir = new File(tmpDir, "out");
      String err = checkHandleFlagsErr("-o", outDir.getPath(), "-t", outDir.getPath());
      TestUtils.containsAllUnwrapped(err, "Error: The specified SDF, \"" + outDir.getPath() + "\", does not exist.");

      err = checkHandleFlagsErr("-o", outDir.getPath(), "-t", template.getPath(), "--father", "blah", in.getPath());
      TestUtils.containsAllUnwrapped(err, "Must set --father and --mother flags and at least one of --son or --daughter flags");

      err = checkHandleFlagsErr("-o", outDir.getPath(), "-t", template.getPath(), "--mother", "blah", in.getPath());
      TestUtils.containsAllUnwrapped(err, "Must set --father and --mother flags and at least one of --son or --daughter flags");

      err = checkHandleFlagsErr("-o", outDir.getPath(), "-t", template.getPath(), "--son", "blah", in.getPath());
      TestUtils.containsAllUnwrapped(err, "Must set --father and --mother flags and at least one of --son or --daughter flags");

      err = checkHandleFlagsErr("-o", outDir.getPath(), "-t", template.getPath(), "--daughter", "blah", in.getPath());
      TestUtils.containsAllUnwrapped(err, "Must set --father and --mother flags and at least one of --son or --daughter flags");

      err = checkHandleFlagsErr("-o", outDir.getPath(), "-t", template.getPath(), "--father", "blah", in.getPath(), "-p", tmpFile.getPath());
      TestUtils.containsAllUnwrapped(err, "Cannot use --pedigree in conjunction with --father, --mother, --son, and --daughter flags");

      err = checkHandleFlagsErr("-o", outDir.getPath(), "-t", template.getPath(), "--mother", "blah", in.getPath(), "-p", tmpFile.getPath());
      TestUtils.containsAllUnwrapped(err, "Cannot use --pedigree in conjunction with --father, --mother, --son, and --daughter flags");

      err = checkHandleFlagsErr("-o", outDir.getPath(), "-t", template.getPath(), "--son", "blah", in.getPath(), "-p", tmpFile.getPath());
      TestUtils.containsAllUnwrapped(err, "Cannot use --pedigree in conjunction with --father, --mother, --son, and --daughter flags");

      err = checkHandleFlagsErr("-o", outDir.getPath(), "-t", template.getPath(), "--daughter", "blah", in.getPath(), "-p", tmpFile.getPath());
      TestUtils.containsAllUnwrapped(err, "Cannot use --pedigree in conjunction with --father, --mother, --son, and --daughter flags");

      err = checkHandleFlagsErr("-o", outDir.getPath(), "-t", template.getPath(), "--father", "blah", "--mother", "blah", "--son", "foo", in.getPath());
      TestUtils.containsAllUnwrapped(err, "Father and mother must be different samples");

      err = checkHandleFlagsErr("-o", outDir.getPath(), "-t", template.getPath(), "--father", "blah", "--mother", "foo", "--son", "blah", in.getPath());
      TestUtils.containsAllUnwrapped(err, "Son must be different sample to mother and father");

      err = checkHandleFlagsErr("-o", outDir.getPath(), "-t", template.getPath(), "--father", "foo", "--mother", "blah", "--son", "blah", in.getPath());
      TestUtils.containsAllUnwrapped(err, "Son must be different sample to mother and father");

      err = checkHandleFlagsErr("-o", outDir.getPath(), "-t", template.getPath(), "--father", "blah", "--mother", "foo", "--daughter", "blah", in.getPath());
      TestUtils.containsAllUnwrapped(err, "Daughter must be different sample to mother and father");

      err = checkHandleFlagsErr("-o", outDir.getPath(), "-t", template.getPath(), "--father", "foo", "--mother", "blah", "--daughter", "blah", in.getPath());
      TestUtils.containsAllUnwrapped(err, "Daughter must be different sample to mother and father");

      err = checkHandleFlagsErr("-o", outDir.getPath(), "-t", template.getPath(), "--father", "foo", "--mother", "bar", "--son", "blah", "--son", "blah", in.getPath());
      TestUtils.containsAllUnwrapped(err, "Individual sons must be different sample to other sons");

      err = checkHandleFlagsErr("-o", outDir.getPath(), "-t", template.getPath(), "--father", "foo", "--mother", "bar", "--daughter", "blah", "--daughter", "blah", in.getPath());
      TestUtils.containsAllUnwrapped(err, "Individual daughters must be different sample to other daughters");

      err = checkHandleFlagsErr("-o", outDir.getPath(), "-t", template.getPath(), "--father", "foo", "--mother", "bar", "--son", "blah", "--daughter", "blah", in.getPath());
      TestUtils.containsAllUnwrapped(err, "Son and daughter samples must be different");

      checkHandleFlagsOut("-o", outDir.getPath(), "-t", template.getPath(), "--father", "foo", "--mother", "bar", "--son", "doe", "--daughter", "ray", in.getPath());

      final File in2 = new File(tmpDir, "alignments2.sam.gz");
      FileHelper.stringToGzFile(SAMHEADER + SAM_BODY, in2);
      new TabixIndexer(in2, new File(tmpDir, "alignments2.sam.gz.tbi")).saveSamIndex();

      err = checkMainInitBadFlags("-o", outDir.getPath(), "-t", template.getPath(), "-p", tmpFile.getPath(), in2.getPath());
      TestUtils.containsAllUnwrapped(err, "should contain exactly 6");

      err = checkMainInitWarn("-o", outDir.getPath(), "-t", template.getPath(), "-p", tmpFile2.getPath(), in.getPath(), in.getPath(), "--" + AbstractMultisampleCli.NO_CALIBRATION);
      assertTrue(err, err.contains("assuming autosomal inheritance"));
    } finally {
      FileHelper.deleteAll(tmpDir);
      FileHelper.deleteAll(tmpFile);
      FileHelper.deleteAll(tmpFile2);
    }
  }

  public void testLowCoverage() throws Exception {
    try (TestDirectory dir = new TestDirectory()) {
      final File tmpFile2 = FileUtils.stringToFile(FAMILY_PED, new File(dir, "pedfile"));
      final File in = new File(dir, "alignments.sam.gz");
      FileHelper.stringToGzFile(SAM_FAMILY, in);
      new TabixIndexer(in, new File(dir, "alignments.sam.gz.tbi")).saveSamIndex();

      final File template = new File(dir, "template");
      ReaderTestUtils.getDNADir(">g1\nacgtacgtacgtacgtacgt", template);
      ReferenceTextBuilder.createDiploid().addSequence("g1", Sex.EITHER, Ploidy.DIPLOID, true).writeToSdfDir(template);

      final MainResult calRes = MainResult.run(new RecalibrateCli(), "-t", template.getPath(), in.getPath());
      assertEquals(calRes.err(), 0, calRes.rc());

      final File outDir = new File(dir, "out");
      final String err = checkMainInitWarn("-o", outDir.getPath(), "-t", template.getPath(), "-p", tmpFile2.getPath(), in.getPath(), in.getPath());
      TestUtils.containsAll(err, "Calibration indicates very low coverage", "Check that appropriate target regions were supplied during mapping/calibration");
    }
  }

  public void testPruneHypothesis() throws Exception {
    final File tmpDir = FileHelper.createTempDirectory();
    final File tmpFile = FileUtils.stringToFile(FAMILY_PED, FileHelper.createTempFile());
    try {
      final File in = new File(tmpDir, "alignments.sam.gz");
      FileHelper.stringToGzFile(SAM_FAMILY, in);
      new TabixIndexer(in, new File(tmpDir, "alignments.sam.gz.tbi")).saveSamIndex();
      final File template = new File(tmpDir, "template");
      ReaderTestUtils.getReaderDNA(">g1\nacgtacgtacgtacgtacgt", template, null);
      final File outDir = new File(tmpDir, "out");
      final File in2 = new File(tmpDir, "alignments2.sam.gz");
      FileHelper.stringToGzFile(SAMHEADER + SAM_BODY, in2);
      new TabixIndexer(in2, new File(tmpDir, "alignments2.sam.gz.tbi")).saveSamIndex();

      final VariantParams p = checkMakeParamsOut("-o", outDir.getPath(), "-t", template.getPath(), "-p", tmpFile.getPath(), in.getPath(), in.getPath(), "--" + AbstractMultisampleCli.NO_CALIBRATION);
      assertEquals(true, p.pruneHypotheses());
    } finally {
      FileHelper.deleteAll(tmpDir);
      FileHelper.deleteAll(tmpFile);
    }
  }
}
