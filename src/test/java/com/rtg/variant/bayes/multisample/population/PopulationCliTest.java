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
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;
import com.rtg.variant.VariantParams;
import com.rtg.variant.bayes.multisample.AbstractMultisampleCli;
import com.rtg.variant.bayes.multisample.family.FamilyCliTest;

/**
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
    try (final TestDirectory tmpDir = new TestDirectory()) {
      final File tmpFile = FileUtils.stringToFile("original-derived TEST cancer contamination=0.13", FileHelper.createTempFile(tmpDir));
      final File tmpFile2 = FileUtils.stringToFile("0\tfather\t0\t0\t1\t0\n0\tmother\t0\t0\t2\t0\n0\tchild\tfather\tmother\t1\t0\n", FileHelper.createTempFile(tmpDir));
      checkValidator(tmpDir, tmpFile, tmpFile2, SharedSamConstants.SAM_FAMILY);
      checkReferenceWarning(tmpDir, tmpFile2);
    }
  }

  private void checkReferenceWarning(File tmpDir, File tmpFile2) throws Exception {
      final File in = new File(tmpDir, "alignments.sam.gz");
      FileHelper.stringToGzFile(SharedSamConstants.SAM_FAMILY, in);
      new TabixIndexer(in, new File(tmpDir, "alignments.sam.gz.tbi")).saveSamIndex();

      final File template = new File(tmpDir, "template");
      ReaderTestUtils.getDNADir(">g1\nacgtacgtacgtacgtacgt", template);
      final File outDir = new File(tmpDir, "out");

      final String err = checkMainInitWarn("-o", outDir.getPath(), "-t", template.getPath(), "--pedigree", tmpFile2.getPath(), in.getPath(), in.getPath(), "--" + AbstractMultisampleCli.NO_CALIBRATION);
      assertTrue(err, err.contains("assuming autosomal inheritance"));
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

    checkHandleFlagsOut("-o", outDir.getPath(), "-t", template.getPath(), "--pedigree", tmpFile2.getPath(), in.getPath(), in.getPath(), "--" + AbstractMultisampleCli.NO_CALIBRATION);
  }

  @Override
  protected ParamsCli<VariantParams> getParamsCli() {
    return new PopulationCli();
  }
}
