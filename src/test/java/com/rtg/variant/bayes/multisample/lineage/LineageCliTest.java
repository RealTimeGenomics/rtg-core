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
package com.rtg.variant.bayes.multisample.lineage;

import static com.rtg.util.StringUtils.LS;

import java.io.ByteArrayInputStream;
import java.io.File;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.MainResult;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.sam.SharedSamConstants;
import com.rtg.tabix.TabixIndexer;
import com.rtg.util.StringUtils;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.BgzipFileHelper;
import com.rtg.util.test.FileHelper;
import com.rtg.variant.bayes.multisample.AbstractMultisampleCli;
import com.rtg.variant.bayes.multisample.population.PopulationCliTest;

/**
 */
public class LineageCliTest extends PopulationCliTest {
  // Can't extend AbstractCallerCliTest because it doesn't handle additional required flags

  private static final String EXP_F1 = "Error: You must provide values for -o DIR -p FILE -t SDF" + LS;
  private static final String RESOURCES_DIR = "com/rtg/variant/bayes/multisample/lineage/resources/";

  @Override
  public String getExpectedF1() {
    return EXP_F1;
  }

  @Override
  public void testInitParams() {
    checkHelp("lineage [OPTION]... -o DIR -p FILE -t SDF FILE+",
        "[OPTION]... -o DIR -p FILE -t SDF -I FILE",
        "Performs a combined cell lineage variant analysis."
    );
  }


  @Override
  public void testValidator() throws Exception {
    final File tmpDir = FileHelper.createTempDirectory();
    final File tmpFile = FileUtils.stringToFile("original-derived TEST cancer contamination=0.13", FileHelper.createTempFile());
    final File tmpFile2 = FileUtils.stringToFile("0\toriginal\t0\t0\t1\t0\n0\tleft\toriginal\t0\t1\t0\n0\tright\toriginal\t0\t1\t0\n", FileHelper.createTempFile());
    try {
      checkValidator(tmpDir, tmpFile, tmpFile2, SharedSamConstants.SAM_LINEAGE);
    } finally {
      FileHelper.deleteAll(tmpDir);
      FileHelper.deleteAll(tmpFile);
      FileHelper.deleteAll(tmpFile2);
    }
  }
  public void testNanoTwoSamples() throws Exception {
    check("ref.fasta", "twoSample.sam", "twoSample.ped", "twoSample_expected.txt");
  }

  public void testNanoTwoSamplesDisconnected() throws Exception {
    check("ref.fasta", "twoSampleDisconnected.sam", "twoSampleDisconnected.ped", "twoSampleDisconnected_expected.txt");
  }

  public void testNanoThreeSamples() throws Exception {
    check("ref.fasta", "threeSample.sam", "threeSample.ped", "threeSample_expected.txt");
  }

  void check(String ref, String sam, String ped, String expected) throws Exception {
    try (TestDirectory tmp = new TestDirectory()) {
      final File reference = new File(tmp, "template");
      final String refString = FileHelper.resourceToString(RESOURCES_DIR + ref);
      ReaderTestUtils.getDNADir(refString, reference);
      final String samContents = FileHelper.resourceToString(RESOURCES_DIR + sam);
      final File pedigree = new File(tmp, "pedigree.ped");
      FileHelper.resourceToFile(RESOURCES_DIR + ped, pedigree);
      final File samFile = createSam(samContents, tmp, "lineage");
      final File output = new File(tmp, "output");

      final MainResult res = MainResult.run(new LineageCli(), "--output", output.toString()
        , "--pedigree", pedigree.toString()
        , "--template", reference.toString()
        , samFile.toString()
        , "--" + AbstractMultisampleCli.NO_CALIBRATION);
      assertEquals(res.err(), 0, res.rc());
      final String result = StringUtils.grepMinusV(FileHelper.zipFileToString(new File(output, "snps.vcf.gz")), "^#");
      mNano.check(expected, result);
    }
  }

  private File createSam(final String sam, final File dir, final String name) throws Exception {
    final File file = new File(dir, name + ".sam.gz");
    BgzipFileHelper.streamToBgzipFile(new ByteArrayInputStream(sam.getBytes()), file);
    new TabixIndexer(file, new File(dir, file.getName() + ".tbi")).saveSamIndex();
    return file;
  }

  @Override
  protected AbstractCli getCli() {
    return new LineageCli();
  }
}
