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

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.util.Properties;

import com.rtg.launcher.AbstractNanoTest;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.tabix.TabixIndexer;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.PropertiesUtils;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.Utils;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.BgzipFileHelper;
import com.rtg.util.test.FileHelper;
import com.rtg.util.test.NanoRegression;
import com.rtg.variant.bayes.multisample.AbstractMultisampleCli;

/**
 */
public class PopulationNanoTest extends AbstractNanoTest {

  private static final String RESOURCES_DIR = "com/rtg/variant/bayes/multisample/family/resources/";
  /** the number of expected bases in denovo tests */
  public static final int FAMILY_DENOVO_METRIC = 9400;

  /** Pedigree file text */
  public static final String FAMILY_PED = ""
      + "0\tsm_dad\t0\t0\t1\t0" + StringUtils.LS
      + "0\tsm_mom\t0\t0\t2\t0" + StringUtils.LS
      + "0\tsm_son1\tsm_dad\tsm_mom\t1\t0" + StringUtils.LS
      + "0\tsm_son2\tsm_dad\tsm_mom\t1\t0" + StringUtils.LS
      + "0\tsm_daughter1\tsm_dad\tsm_mom\t2\t0" + StringUtils.LS
      + "0\tsm_daughter2\tsm_dad\tsm_mom\t2\t0" + StringUtils.LS;

  /** Reference sequence */
  public static final String DENOVO_REF = ">simulatedSequence1 100 1" + StringUtils.LS
      + "CTTGTATCGTGTGTAATGAAAGCGCACTCATCATAGTTCTTTTGAGTCGATACTTAGCTCCCATGTGATGCGGATGACCGGTGAACC"
      + "GCAATATACATGAGTTTTAGGAACGGTACTAGTACCTCTAATACGCACTACTTTGTTCGATTACCATGGATAATCCGCCCTTGTTGG"
      + "ATTGCAATCGAACGGCAGGATCGCT" + StringUtils.LS;

  public void testNoDenovo() throws Exception {
    final String sam = FileHelper.resourceToString(RESOURCES_DIR + "family_denovo.sam");
    try (TestDirectory dir = new TestDirectory()) {
      final File priors = makeDenovoPriors(dir, 0.0, 0.0);
      check(mNano, DENOVO_REF, sam, FAMILY_PED, "nodenovo", new String[] {"--Xpriors", priors.getPath(), "--Xcontrary-probability", "1", "--Xmax-em-iterations", "50"}, FAMILY_DENOVO_METRIC, false);
    }
  }
  public static File makeDenovoPriors(File dir, double reference, double nonReference) throws IOException, InvalidParamsException {
    final Properties p = PropertiesUtils.getPriorsResource("human", PropertiesUtils.PropertyType.PRIOR_PROPERTY);
    p.setProperty("denovo_reference_rate", Double.toString(reference));
    p.setProperty("denovo_non_reference_rate", Double.toString(nonReference));
    final File testDenovoPriors = new File(dir, "testDenovoPriors");
    try (Writer writer = new FileWriter(testDenovoPriors)) {
      p.store(writer, "PopulationNanoTest properties");
    }
    return testDenovoPriors;
  }

  public void testDenovo() throws Exception {
    final String sam = FileHelper.resourceToString(RESOURCES_DIR + "family_denovo.sam");
    try (TestDirectory dir = new TestDirectory()) {
      final File priors = makeDenovoPriors(dir, 0.002, 0.002);
      check(mNano, DENOVO_REF, sam, FAMILY_PED, "denovo", new String[] {"--Xpriors", priors.getPath(), "--Xcontrary-probability", "1", "--Xmax-em-iterations", "50"}, FAMILY_DENOVO_METRIC, false);
    }
  }

  public void testDenovoHomozygousOnly() throws Exception {
    final String sam = FileHelper.resourceToString(RESOURCES_DIR + "family_denovo.sam");
    try (TestDirectory dir = new TestDirectory()) {
      final File priors = makeDenovoPriors(dir, 0.002, 0);
      check(mNano, DENOVO_REF, sam, FAMILY_PED, "denovo_hom", new String[] {"--Xpriors", priors.getPath(), "--Xcontrary-probability", "1", "--Xmax-em-iterations", "50"}, FAMILY_DENOVO_METRIC, false);
    }
  }

  public void testDenovoMixed() throws Exception {
    final String sam = FileHelper.resourceToString(RESOURCES_DIR + "family_denovo.sam");
    try (TestDirectory dir = new TestDirectory()) {
      final File priors = makeDenovoPriors(dir, 0.002, 0.000002);
      check(mNano, DENOVO_REF, sam, FAMILY_PED, "denovo_mixed", new String[] {"--Xpriors", priors.getPath(), "--Xcontrary-probability", "1", "--Xmax-em-iterations", "50"}, FAMILY_DENOVO_METRIC, false);
    }
  }

  public void testProp() throws Exception {
    final String sam = FileHelper.resourceToString(RESOURCES_DIR + "family_denovo.sam");
    try (TestDirectory dir = new TestDirectory()) {
      final File priors = makeDenovoPriors(dir, 0.0, 0.0);
      check(mNano, DENOVO_REF, sam, FAMILY_PED, "nodenovo_prop", new String[] {"--Xpriors", priors.getPath(), "--Xprop-priors", "--Xmax-em-iterations", "50"}, FAMILY_DENOVO_METRIC, false);
    }
    try (TestDirectory dir = new TestDirectory()) {
      final File priors = makeDenovoPriors(dir, 0.002, 0.002);
      check(mNano, DENOVO_REF, sam, FAMILY_PED, "denovo_prop", new String[] {"--Xpriors", priors.getPath(), "--Xprop-priors", "--Xmax-em-iterations", "50"}, FAMILY_DENOVO_METRIC, false);
    }
  }


  /**
   * Evaluate a {@code NanoRegression} which uses the population caller.
   * public so it can be used from multi sample task to test indel triggering
   * @param template template sequence
   * @param samString sam file content
   * @param ped pedigree file
   * @param expectedPrefix regression name
   * @param extraArgs arguments to add to command line
   * @param metric number of records expected to be processed
   * @param bed test an expected bed file
   * @throws Exception when IO explodes
   */
  public static void check(NanoRegression nano, final String template, final String samString, final String ped, final String expectedPrefix, String[] extraArgs, int metric, boolean bed) throws Exception {
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
          "--" + AbstractMultisampleCli.NO_CALIBRATION
      };

      final String[] args;
      if (extraArgs != null) {
        args = Utils.append(baseArgs, extraArgs);
      } else {
        args = baseArgs;
      }

      final MemoryPrintStream ps = new MemoryPrintStream();

      final PopulationCli cli = new PopulationCli();
      final int code = cli.mainInit(args, new ByteArrayOutputStream(), ps.printStream());
      assertEquals(ps.toString(), 0, code);

      final String result = FileUtils.fileToString(new File(out, "snps.vcf"));
      final String actualFixed = TestUtils.stripVcfHeader(result);

      //System.err.println(actualFixed);
      //final String expectedFinal = FileHelper.resourceToString("com/rtg/variant/bayes/multisample/cancer/resources/somaticnanotest" + expectedPrefix + ".vcf").replaceAll("\r", "");
      //assertEquals(expectedFinal, actualFixed.replaceAll("\r", ""));
      nano.check("populationnanotest-" + expectedPrefix + ".vcf", actualFixed, false);
      final String usageLog = cli.usageLog();
      //System.err.println(usageLog);
      TestUtils.containsAll(usageLog, "[Usage beginning module=population runId=", ", Usage end module=population runId=", " metric=" + metric + " success=true]");
      if (bed) {
        final String actualBedFixed = FileUtils.fileToString(new File(out, "regions.bed"));
        nano.check("populationnanotest-" + expectedPrefix + ".bed", actualBedFixed, false);
      }
    } finally {
      assertTrue(FileHelper.deleteAll(tmp));
    }
  }

  private static File createIndexedSamFile(final String cancer, final File tmp, final String tag) throws Exception {
    final File file = new File(tmp, tag + ".sam.gz");
    BgzipFileHelper.streamToBgzipFile(new ByteArrayInputStream(cancer.getBytes()), file);
    new TabixIndexer(file, new File(tmp, tag + ".sam.gz.tbi")).saveSamIndex();
    return file;
  }
}
