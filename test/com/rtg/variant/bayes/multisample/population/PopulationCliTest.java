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

package com.rtg.variant.bayes.multisample.population;

import static com.rtg.util.StringUtils.LS;

import java.io.File;
import java.io.IOException;

import com.rtg.launcher.AbstractParamsCliTest;
import com.rtg.launcher.ParamsCli;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.sam.SharedSamConstants;
import com.rtg.tabix.TabixIndexer;
import com.rtg.tabix.UnindexableDataException;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;
import com.rtg.variant.VariantParams;
import com.rtg.variant.bayes.multisample.AbstractMultisampleCli;
import com.rtg.variant.bayes.multisample.family.FamilyCliTest;

/**
 *         Date: 13/03/12
 *         Time: 10:28 AM
 */
public class PopulationCliTest extends AbstractParamsCliTest<VariantParams> {

  private static final String EXP_F1 = "Error: You must provide values for -o DIR -p FILE -t SDF" + LS;
  public String getExpectedF1() {
    return EXP_F1;
  }

  public void testErrorF1() {
    final String err = checkHandleFlagsErr();
    final String expF1 = getExpectedF1();
    assertTrue("<" + expF1 + "> was not contained in <" + err + ">", err.contains(expF1));
  }

  public void testInitParams() {
    checkHelp("population [OPTION]... -o DIR -p FILE -t SDF FILE+",
        "[OPTION]... -o DIR -p FILE -t SDF -I FILE",
        "Performs a multiple sample variant analysis of many, potentially related, genomes"
    );
  }

  public void testValidator() throws Exception {
    final File tmpDir = FileHelper.createTempDirectory();
    final File tmpFile = FileUtils.stringToFile("original-derived TEST cancer contamination=0.13", FileHelper.createTempFile());
    final File tmpFile2 = FileUtils.stringToFile("0\tfather\t0\t0\t1\t0\n0\tmother\t0\t0\t2\t0\n0\tchild\tfather\tmother\t1\t0\n", FileHelper.createTempFile());
    try {
      checkValidator(tmpDir, tmpFile, tmpFile2, SharedSamConstants.SAM_FAMILY);

    } finally {
      FileHelper.deleteAll(tmpDir);
      FileHelper.deleteAll(tmpFile);
      FileHelper.deleteAll(tmpFile2);
    }
  }

  public void checkValidator(File tmpDir, File tmpFile, File tmpFile2, String sam) throws IOException, UnindexableDataException {
    final File in = new File(tmpDir, "alignments.sam.gz");
    FileHelper.stringToGzFile(sam, in);
    new TabixIndexer(in, new File(tmpDir, "alignments.sam.gz.tbi")).saveSamIndex();

    final File template = new File(tmpDir, "template");
    ReaderTestUtils.getDNADir(">g1\nacgtacgtacgtacgtacgt", template);
    final File outDir = new File(tmpDir, "out");

    String err = checkHandleFlagsErr("-o", outDir.getPath(), "-t", outDir.getPath(), "--pedigree", tmpFile.getPath());
    assertTrue(err, err.contains("Error: The specified SDF, \"" + outDir.getPath() + "\", does not exist."));

    err = checkHandleFlagsErr("-o", outDir.getPath(), "-t", template.getPath(), tmpFile.getPath());
    assertTrue(err, err.contains("You must provide a value for -p FILE"));

    err = checkHandleFlagsErr("-o", outDir.getPath(), "-t", template.getPath(), "--pedigree", tmpFile.getPath());
    assertTrue(err, err.contains("No input files specified"));

    checkHandleFlagsOut("-o", outDir.getPath(), "-t", template.getPath(), "--pedigree", tmpFile.getPath(), in.getPath());

    final File in2 = new File(tmpDir, "alignments2.sam.gz");
    FileHelper.stringToGzFile(FamilyCliTest.SAMHEADER + SharedSamConstants.SAM_BODY, in2);
    new TabixIndexer(in2, new File(tmpDir, "alignments2.sam.gz.tbi")).saveSamIndex();

    err = checkMainInitBadFlags("-o", outDir.getPath(), "-t", template.getPath(), "--pedigree", tmpFile.getPath(), in2.getPath());
    assertTrue(err, err.contains("should contain exactly 6 columns"));

    err = checkHandleFlagsErr("-o", outDir.getPath(), "-t", template.getPath(), "--pedigree", tmpFile.getPath(), in.getPath(), "--min-avr-score", "1.1");
    assertTrue(err, err.contains("--min-avr-score must be in the range [0.0, 1.0]"));

    err = checkHandleFlagsErr("-o", outDir.getPath(), "-t", template.getPath(), "--pedigree", tmpFile.getPath(), in.getPath(), "--min-avr-score", "-0.1");
    assertTrue(err, err.contains("--min-avr-score must be in the range [0.0, 1.0]"));

    err = checkMainInitWarn("-o", outDir.getPath(), "-t", template.getPath(), "--pedigree", tmpFile2.getPath(), in.getPath(), in.getPath(), "--" + AbstractMultisampleCli.NO_CALIBRATION);
    assertTrue(err, err.contains("assuming autosomal inheritance"));
  }

  @Override
  protected ParamsCli<VariantParams> getParamsCli() {
    return new PopulationCli();
  }
}
