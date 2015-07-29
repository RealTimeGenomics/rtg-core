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

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;

import com.rtg.reader.ReaderTestUtils;
import com.rtg.tabix.TabixIndexer;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.Talkback;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.BgzipFileHelper;
import com.rtg.util.test.FileHelper;
import com.rtg.util.test.NanoRegression;
import com.rtg.variant.bayes.multisample.population.PopulationNanoTest;

import junit.framework.TestCase;

/**
 */
public class FamilyNanoTest extends TestCase {

  private static final String RESOURCES_DIR = "com/rtg/variant/bayes/multisample/family/resources/";
  private NanoRegression mNano = null;

  @Override
  public void setUp() {
    Diagnostic.setLogStream();
    mNano = new NanoRegression(this.getClass(), false);
  }

  @Override
  public void tearDown() throws Exception {
    // clear the module name so later tests don't report SlimException to the Talkback system
    Talkback.setModuleName(null);
    try {
      mNano.finish();
    } finally {
      mNano = null;
    }
  }

  private static final String FAMILY_PED = ""
      + "0\tsm_dad\t0\t0\t1\t0" + StringUtils.LS
      + "0\tsm_mom\t0\t0\t2\t0" + StringUtils.LS
      + "0\tsm_son1\tsm_dad\tsm_mom\t1\t0" + StringUtils.LS
      + "0\tsm_son2\tsm_dad\tsm_mom\t1\t0" + StringUtils.LS
      + "0\tsm_daughter1\tsm_dad\tsm_mom\t2\t0" + StringUtils.LS
      + "0\tsm_daughter2\tsm_dad\tsm_mom\t2\t0" + StringUtils.LS;

  private static final String DENOVO_REF = ">simulatedSequence1 100 1" + StringUtils.LS
      + "CTTGTATCGTGTGTAATGAAAGCGCACTCATCATAGTTCTTTTGAGTCGATACTTAGCTCCCATGTGATGCGGATGACCGGTGAACC"
      + "GCAATATACATGAGTTTTAGGAACGGTACTAGTACCTCTAATACGCACTACTTTGTTCGATTACCATGGATAATCCGCCCTTGTTGG"
      + "ATTGCAATCGAACGGCAGGATCGCT" + StringUtils.LS;

  public void testNoDenovo() throws Exception {
    final String sam = FileHelper.resourceToString(RESOURCES_DIR + "family_denovo.sam");

    try (TestDirectory dir = new TestDirectory()) {
      final File priors = PopulationNanoTest.makeDenovoPriors(dir, 0.0, 0.0);
      check(DENOVO_REF, sam, FAMILY_PED, "nodenovo", new String[] {"--Xpriors", priors.getPath()});
    }
  }

  public void testDenovo() throws Exception {
    final String sam = FileHelper.resourceToString(RESOURCES_DIR + "family_denovo.sam");
    try (TestDirectory dir = new TestDirectory()) {
      final File priors = PopulationNanoTest.makeDenovoPriors(dir, 0.002, 0.002);
      check(DENOVO_REF, sam, FAMILY_PED, "denovo", new String[] {"--Xpriors", priors.getPath()});
    }
  }

  public void testDenovoHom() throws Exception {
    final String sam = FileHelper.resourceToString(RESOURCES_DIR + "family_denovo.sam");
    try (TestDirectory dir = new TestDirectory()) {
      final File priors = PopulationNanoTest.makeDenovoPriors(dir, 0.002, 0.0);
      check(DENOVO_REF, sam, FAMILY_PED, "denovo_hom", new String[] {"--Xpriors", priors.getPath()});
    }
  }

  private void check(final String template, final String samString, final String ped, final String expectedPrefix, String[] extraArgs) throws Exception {
    //Diagnostic.setLogStream(System.err);
    final File tmp = FileUtils.createTempDir("familynano", "tests");
    try {
      final File templateDir = new File(tmp, "template");
      ReaderTestUtils.getDNADir(template, templateDir);
      final File relationsFile = new File(tmp, "family.ped");
      final File out = new File(tmp, "output");
      FileUtils.stringToFile(ped, relationsFile);
      final File sam = createIndexedSamFile(samString, tmp, "total");


      final String[] baseArgs = {
          "-t", templateDir.getPath(),
          "-p", relationsFile.getPath(),
          "-o", out.getPath(),
          "-Z",
          sam.getPath(),
      };

      final String[] args;
      if (extraArgs != null) {
        args = TestUtils.append(baseArgs, extraArgs);
      } else {
        args = baseArgs;
      }

      final MemoryPrintStream ps = new MemoryPrintStream();

      final FamilyCli cli = new FamilyCli();
      final int code = cli.mainInit(args, new ByteArrayOutputStream(), ps.printStream());
      assertEquals(ps.toString(), 0, code);

      final String result = FileUtils.fileToString(new File(out, "snps.vcf"));
      final String actualFixed = TestUtils.stripVcfHeader(result);

      //System.err.println(actualFixed);
      //final String expectedFinal = FileHelper.resourceToString("com/rtg/variant/bayes/multisample/cancer/resources/somaticnanotest" + expectedPrefix + ".vcf").replaceAll("\r", "");
      //assertEquals(expectedFinal, actualFixed.replaceAll("\r", ""));
      mNano.check("familynanotest-" + expectedPrefix + ".vcf", actualFixed, false);
      final String usageLog = cli.usageLog();
      //System.err.println(usageLog);
      TestUtils.containsAll(usageLog, "[Usage beginning module=family runId=", ", Usage end module=family runId=", " metric=9400 success=true]");
    } finally {
      assertTrue(FileHelper.deleteAll(tmp));
    }
  }

  private File createIndexedSamFile(final String cancer, final File tmp, final String tag) throws Exception {
    final File file = new File(tmp, tag + ".sam.gz");
    BgzipFileHelper.streamToBgzipFile(new ByteArrayInputStream(cancer.getBytes()), file);
    new TabixIndexer(file, new File(tmp, tag + ".sam.gz.tbi")).saveSamIndex();
    return file;
  }

  //To replicate bug 1550
  public void testBug1550() throws IOException {
    //while (true) {
    try (final TestDirectory dir = new TestDirectory()) {
      final File bam = FileHelper.resourceToFile(RESOURCES_DIR + "input.bam", new File(dir, "input.bam"));
      FileHelper.resourceToFile(RESOURCES_DIR + "input.bam.bai", new File(dir, "input.bam.bai"));
      FileHelper.resourceToFile(RESOURCES_DIR + "input.bam.calibration", new File(dir, "input.bam.calibration"));
      final File ped = FileHelper.resourceToFile(RESOURCES_DIR + "WMD.ped", new File(dir, "WMD.ped"));
      final File template = ReaderTestUtils.getDNADir(FileHelper.resourceToString(RESOURCES_DIR + "subset_template.fasta"), new File(dir, "subset_template"));
      final File outputDir = new File(dir, "output");
      final String[] args = {
          "-t", template.getPath(),
          "-o", outputDir.getPath(),
          "-p", ped.getPath(),
          bam.getPath(),
          "--region", "GL000220.1",
          "--Xignore-incompatible-sam-headers" //bams are annoying to clean up, the input contains sequence dictionary entries not in the reference :(
      };
      final FamilyCli cli = new FamilyCli();
      final MemoryPrintStream ps = new MemoryPrintStream();
      final int code = cli.mainInit(args, ps.outputStream(), ps.printStream());
      assertEquals(ps.toString(), 0, code);
    }
    //}
  }
}
