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
package com.rtg.variant.bayes.multisample.family;

import static com.rtg.sam.SharedSamConstants.SAM_BODY;
import static com.rtg.sam.SharedSamConstants.SAM_FAMILY;
import static com.rtg.util.StringUtils.LS;

import java.io.File;

import com.rtg.launcher.AbstractParamsCliTest;
import com.rtg.launcher.ParamsCli;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.tabix.TabixIndexer;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;
import com.rtg.variant.VariantParams;

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
      assertTrue(err, err.contains("Error: The specified SDF, \"" + outDir.getPath() + "\", does not exist."));

      err = checkHandleFlagsErr("-o", outDir.getPath(), "-t", template.getPath(), "--father", "blah", in.getPath());
      assertTrue(err, err.contains("Must set --father and --mother flags and at least one of --son or --daughter flags"));

      err = checkHandleFlagsErr("-o", outDir.getPath(), "-t", template.getPath(), "--mother", "blah", in.getPath());
      assertTrue(err, err.contains("Must set --father and --mother flags and at least one of --son or --daughter flags"));

      err = checkHandleFlagsErr("-o", outDir.getPath(), "-t", template.getPath(), "--son", "blah", in.getPath());
      assertTrue(err, err.contains("Must set --father and --mother flags and at least one of --son or --daughter flags"));

      err = checkHandleFlagsErr("-o", outDir.getPath(), "-t", template.getPath(), "--daughter", "blah", in.getPath());
      assertTrue(err, err.contains("Must set --father and --mother flags and at least one of --son or --daughter flags"));

      err = checkHandleFlagsErr("-o", outDir.getPath(), "-t", template.getPath(), "--father", "blah", in.getPath(), "-p", tmpFile.getPath());
      assertTrue(err, err.contains("Cannot use --pedigree in conjunction with --father, --mother, --son, and --daughter flags"));

      err = checkHandleFlagsErr("-o", outDir.getPath(), "-t", template.getPath(), "--mother", "blah", in.getPath(), "-p", tmpFile.getPath());
      assertTrue(err, err.contains("Cannot use --pedigree in conjunction with --father, --mother, --son, and --daughter flags"));

      err = checkHandleFlagsErr("-o", outDir.getPath(), "-t", template.getPath(), "--son", "blah", in.getPath(), "-p", tmpFile.getPath());
      assertTrue(err, err.contains("Cannot use --pedigree in conjunction with --father, --mother, --son, and --daughter flags"));

      err = checkHandleFlagsErr("-o", outDir.getPath(), "-t", template.getPath(), "--daughter", "blah", in.getPath(), "-p", tmpFile.getPath());
      assertTrue(err, err.contains("Cannot use --pedigree in conjunction with --father, --mother, --son, and --daughter flags"));

      err = checkHandleFlagsErr("-o", outDir.getPath(), "-t", template.getPath(), "--father", "blah", "--mother", "blah", "--son", "foo", in.getPath());
      assertTrue(err, err.contains("Father and mother must be different samples"));

      err = checkHandleFlagsErr("-o", outDir.getPath(), "-t", template.getPath(), "--father", "blah", "--mother", "foo", "--son", "blah", in.getPath());
      assertTrue(err, err.contains("Son must be different sample to mother and father"));

      err = checkHandleFlagsErr("-o", outDir.getPath(), "-t", template.getPath(), "--father", "foo", "--mother", "blah", "--son", "blah", in.getPath());
      assertTrue(err, err.contains("Son must be different sample to mother and father"));

      err = checkHandleFlagsErr("-o", outDir.getPath(), "-t", template.getPath(), "--father", "blah", "--mother", "foo", "--daughter", "blah", in.getPath());
      assertTrue(err, err.contains("Daughter must be different sample to mother and father"));

      err = checkHandleFlagsErr("-o", outDir.getPath(), "-t", template.getPath(), "--father", "foo", "--mother", "blah", "--daughter", "blah", in.getPath());
      assertTrue(err, err.contains("Daughter must be different sample to mother and father"));

      err = checkHandleFlagsErr("-o", outDir.getPath(), "-t", template.getPath(), "--father", "foo", "--mother", "bar", "--son", "blah", "--son", "blah", in.getPath());
      assertTrue(err, err.contains("Individual sons must be different sample to other sons"));

      err = checkHandleFlagsErr("-o", outDir.getPath(), "-t", template.getPath(), "--father", "foo", "--mother", "bar", "--daughter", "blah", "--daughter", "blah", in.getPath());
      assertTrue(err, err.contains("Individual daughters must be different sample to other daughters"));

      err = checkHandleFlagsErr("-o", outDir.getPath(), "-t", template.getPath(), "--father", "foo", "--mother", "bar", "--son", "blah", "--daughter", "blah", in.getPath());
      assertTrue(err, err.contains("Son and daughter samples must be different"));

      checkHandleFlagsOut("-o", outDir.getPath(), "-t", template.getPath(), "--father", "foo", "--mother", "bar", "--son", "doe", "--daughter", "ray", in.getPath());

      final File in2 = new File(tmpDir, "alignments2.sam.gz");
      FileHelper.stringToGzFile(SAMHEADER + SAM_BODY, in2);
      new TabixIndexer(in2, new File(tmpDir, "alignments2.sam.gz.tbi")).saveSamIndex();

      err = checkMainInitBadFlags("-o", outDir.getPath(), "-t", template.getPath(), "-p", tmpFile.getPath(), in2.getPath());
      assertTrue(err, err.contains("should contain exactly 6"));

      checkMainInitOk("-o", outDir.getPath(), "-t", template.getPath(), "-p", tmpFile2.getPath(), in.getPath(), in.getPath());
    } finally {
      FileHelper.deleteAll(tmpDir);
      FileHelper.deleteAll(tmpFile);
      FileHelper.deleteAll(tmpFile2);
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

      final VariantParams p = checkMakeParamsOut("-o", outDir.getPath(), "-t", template.getPath(), "-p", tmpFile.getPath(), in.getPath(), in.getPath());
      assertEquals(true, p.pruneHypotheses());
    } finally {
      FileHelper.deleteAll(tmpDir);
      FileHelper.deleteAll(tmpFile);
    }
  }
}
