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

package com.rtg.variant.bayes.multisample.cancer;


import static com.rtg.util.StringUtils.LS;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.File;

import com.rtg.reader.ReaderTestUtils;
import com.rtg.tabix.TabixIndexer;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.Talkback;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.test.BgzipFileHelper;
import com.rtg.util.test.FileHelper;
import com.rtg.util.test.NanoRegression;

import junit.framework.TestCase;

/**
 */
public class SomaticNanoTest extends TestCase {

  private static final String RESOURCES_DIR = "com/rtg/variant/bayes/multisample/cancer/resources/";
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

  private static final String REF_TEST1 = ">simulatedSequence1" + LS
    + "GTGGAAGAGCCTCTCCTAGGATACCATGCGGCTAGGCTAGGCTCGCCCGGCGGTTCCAGGAGCTGGCGAGTGCC"
    + "TCGTTTCCTTACGTGGGACGCTCAAATTCTGCTCTGGTTGTTTAAACTGAATTAGGGGAAGCTGCTAGCGAACT"
    + "CTGCCCAAAACAAGAAACAACCCCGTCCTCTTACGCCGTAACCAGCCA";


  //test an interesting call using cancer caller (homozygous SNP)
  public void test1() throws Exception {
    checkCancer("1", "1", REF_TEST1, 6550L, "--all", "--keep-duplicates", "--Xinclude-germline", "--loh", "0.1");
  }

  public void testGainOfReference() throws Exception {
    checkCancer("6", "6", REF_TEST1, 2680L, "--all", "--keep-duplicates", "--loh", "0.1", "-G");
  }

  private static final String REF_TEST2 = ">test1" + LS
      + "TCGTACGTCT";

  //test an interesting call using cancer caller (homozygous add)
  // Cancer / normal Aviews
  //   TCGTACG_TCT
  //    CGTACGG
  //    CGTACGG
  //     GTACGGT
  //     GTACGGT
  //     GTACGGT
  //      TACGGTC
  //      TACGGTC
  //       ACGGTCT
  //       ACGGTCT
  public void test2() throws Exception {
    checkCancer("2", "2", REF_TEST2, 126L, "--all", "--keep-duplicates", "--Xinclude-germline", "--loh", "0.1");
  }

  private static final String REF_TEST3 = ">test1" + LS
      + "GCTGCTACTCTC";

  //test an interesting call using cancer caller (homozygous insertion)
  // Cancer / normal Aviews
  //  GCTGCTACTCTC
  //   CTGC-ACT
  //   CTGC-ACT
  //    TGC-ACTC
  //    TGC-ACTC
  //     GC-ACTCT
  //     GC-ACTCT
  //      C-ACTCTC
  //      C-ACTCTC
  public void test3() throws Exception {
    checkCancer("3", REF_TEST3, 112L);
  }

  // test an interesting call using cancer caller (homozygous SNP with skipping
  // first 50 Nts)
  public void test4() throws Exception {
    checkCancer("1", "4", REF_TEST1, 2450L, "--all", "--region", "simulatedSequence1:90+10", "--Xinclude-germline", "--loh", "0.1");
  }

  private static final String REF_TEST5 = ""
      + ">test1" + LS
      + "TCGTACGGTCT"
      ;

  // Long add case seen as a bug in wild
  public void test5() throws Exception {
    checkCancer("5", "5", REF_TEST5, 148L, "--keep-duplicates");
  }

  private static final String LONG_INSERT_TEMPLATE = ">chr17" + LS
      + "CATTCCCAAGTCACATGACATCGTTTTGAAACTCTGTCATTCATAAATGGGGCATCCCTTATGACAACAGCGTTATAGTCCGACGGGAGTAAATAAAAACAGCTGTGTTACCATGTCAGCAACCTTGAGGTGCCCATAGGTGAAC";

  // Long add case seen as a bug in wild - this shouldn't pass
  public void testLongInsert() throws Exception {
    checkCancer("LongInsert", "LongInsert", LONG_INSERT_TEMPLATE, 1840L, "--keep-duplicates");
  }

  // Long add case but N removed from normal evidence
  public void testLongInsertNoN() throws Exception {
    checkCancer("LongInsertNoN", "LongInsertNoN", LONG_INSERT_TEMPLATE, 1840L, "--keep-duplicates");
  }

  // Set loh should get the LOH call
  public void testLoh() throws Exception {
    checkCancer("Loh", "Loh", REF_TEST5, 126L, "--keep-duplicates", "--loh", "0.1");
  }

  // Set loh to 0 should not get the LOH call
  public void testNoLoh() throws Exception {
    checkCancer("Loh", "5b", REF_TEST5, 126L, "--keep-duplicates", "--loh", "0.0", "--Xinclude-germline");
  }

  // Set loh should get the LOH call - original is heterozygous
  public void testLohDiploid() throws Exception {
    checkCancer("LohDiploid", "LohDiploid", REF_TEST5, 126L, "--keep-duplicates", "--loh", "0.1");
  }

  // Set loh to 1.0 - an important use case - see testLohDiploid
  public void testLohDiploid1() throws Exception {
    checkCancer("LohDiploid", "LohDiploid1", REF_TEST5, 126L, "--keep-duplicates", "--loh", "1.0");
  }


  private void checkCancer(final String id, final String template, final long usageExp) throws Exception {
    checkCancer(id, id, template, usageExp);
  }

  private void checkCancer(final String idIn, final String idOut, final String template, final long usageExp, final String... options) throws Exception {
    final String relation = "original-derived normal cancer contamination=0.5";
    final String cancer = FileHelper.resourceToString(RESOURCES_DIR + "test" + idIn + "_cancer.sam");
    final String normal = FileHelper.resourceToString(RESOURCES_DIR + "test" + idIn + "_normal.sam");

    check(template, normal, cancer, relation, idOut, options, usageExp);
  }

  public void check(final String template, final String normalStr, final String cancerStr, final String relation, final String expectedPrefix, String[] extraArgs, long usageExp) throws Exception {
    //Diagnostic.setLogStream(System.err);
    final File tmp = FileUtils.createTempDir("somnano", "tests");
    try {
      final File templateDir = new File(tmp, "template");
      ReaderTestUtils.getDNADir(template, templateDir);
      final File relationsFile = new File(tmp, "relations.relations");
      final File out = new File(tmp, "output");
      FileUtils.stringToFile(relation, relationsFile);
      final File cancer = createIndexedSamFile(cancerStr, tmp, "cancer");
      final File normal = createIndexedSamFile(normalStr, tmp, "normal");
      final String[] args1 = {
          "-t", templateDir.getPath(),
          "-r", relationsFile.getPath(),
          "-o", out.getPath(),
          "-Z",
          normal.getPath(),
          cancer.getPath()
      };

      final String[] args;
      if (extraArgs != null) {
        args = TestUtils.append(args1, extraArgs);
      } else {
        args = args1;
      }

      final MemoryPrintStream ps = new MemoryPrintStream();

      final SomaticCli cli = new SomaticCli();
      final int code = cli.mainInit(args, new ByteArrayOutputStream(), ps.printStream());
      assertEquals(ps.toString(), 0, code);

      final String result = FileUtils.fileToString(new File(out, "snps.vcf"));
      final String actual = result.replace("Version", "");
      final String actualFixed = actual.replaceAll("##CL=.*\n", "")
          .replaceAll("##TEMPLATE-SDF-ID=.*\n", "")
          .replaceAll("##RUN-ID=.*\n", "")
          .replaceAll("##fileDate=.*\n", "")
          .replaceAll("##reference=.*\n", "")
          .replaceAll("##INFO=.*\n", "")
          .replaceAll("##SAMPLE=.*\n", "")
          .replaceAll("##FORMAT=.*\n", "")
          .replaceAll("##FILTER=.*\n", "")
          .replaceAll("##PEDIGREE=.*\n", "")
          .replaceAll("##fileformat=.*\n", "")
          .replaceAll("##contig=.*\n", "")
          .replaceAll("##source=.*\n", "");

      //System.err.println(actualFixed);
      //final String expectedFinal = FileHelper.resourceToString("com/rtg/variant/bayes/multisample/cancer/resources/somaticnanotest" + expectedPrefix + ".vcf").replaceAll("\r", "");
      //assertEquals(expectedFinal, actualFixed.replaceAll("\r", ""));
      mNano.check("somaticnanotest" + expectedPrefix + ".vcf", actualFixed, false);
      final String usageLog = cli.usageLog();
      //System.err.println(usageLog);
      TestUtils.containsAll(usageLog, "[Usage beginning module=somatic runId=", ", Usage end module=somatic runId=", " metric=" + usageExp + " success=true]");
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

  // Test situation where N's on ref and evicence is shorter sequence of AAAAA.  See Bug#1591
  public void testIndelWithUnknownRef() throws Exception {
    final String tmpl = ">ref\n"
    + "AGTTAAAGAGTGAAACCCTGATAGTCTTACCCCAAGGCCAAAGTCCTATTTTATTATTTTTATATTCTTACTATATATTATACAAATCTTCATTGCAA"
    + "GTTTNNNNNNNNNNNNNNNNNNNNAAGTAAAAACATAAGAAATCTAATTTTTGTATATAAAAGCTGTAAACTAAATTATATATATACACATACATACA"
    + "TACGTGTGTGTGTGTATATATATATACATATATAACCTATGGATTAGGAAAATTTATTGCTTCAACAAACTAAGGGGATTACTTCCCATAAAATTAGT\n";
    final String sam = FileHelper.resourceToString("com/rtg/variant/resources/n2-mappings.sam");
    final String sam2 = FileHelper.resourceToString("com/rtg/variant/resources/n2-mappings2.sam");
    final String relation = "original-derived BACT_SAMPLE CANCER contamination=0";
    check(tmpl, sam, sam2, relation, "", new String[0], 1386L);
  }

  public void testContaminationBug() throws Exception {
    // Not clear what this is actually testing
    final String tmpl = FileHelper.resourceToString(RESOURCES_DIR + "contaminationbug.ref");
    final String sam = FileHelper.resourceToString(RESOURCES_DIR + "contaminationbug.sam");
    final String relation = "original-derived H318-N H318-T contamination=0";
    check(tmpl, sam, sam, relation, "contaminationbugtest", new String[0], 900L);
  }
}
