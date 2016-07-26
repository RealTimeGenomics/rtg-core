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
package com.rtg.variant.bayes.multisample.lineage;

import static com.rtg.util.StringUtils.LS;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.File;

import com.rtg.launcher.AbstractCli;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.sam.SharedSamConstants;
import com.rtg.tabix.TabixIndexer;
import com.rtg.util.StringUtils;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.MemoryPrintStream;
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

      final LineageCli cli = new LineageCli();
      final String[] args = {
          "--output", output.toString()
          , "--pedigree", pedigree.toString()
          , "--template", reference.toString()
          , samFile.toString()
          , "--" + AbstractMultisampleCli.NO_CALIBRATION
      };
      final MemoryPrintStream ps = new MemoryPrintStream();
      final int code = cli.mainInit(args, new ByteArrayOutputStream(), ps.printStream());
      if (code != 0 || !"".equals(ps.toString())) {
        fail(ps.toString());
      }

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
