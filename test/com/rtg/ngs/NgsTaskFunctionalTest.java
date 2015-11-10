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
package com.rtg.ngs;

import static com.rtg.util.StringUtils.LS;
import static com.rtg.util.StringUtils.TAB;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.io.StringWriter;
import java.util.Collections;

import com.rtg.index.hash.ngs.NgsHashLoopImpl;
import com.rtg.launcher.SequenceParams;
import com.rtg.mode.DnaUtils;
import com.rtg.mode.SequenceMode;
import com.rtg.ngs.NgsTestUtils.OverriddenNgsOutputParams;
import com.rtg.ngs.NgsTestUtils.ParamsParams;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.reader.SequencesReader;
import com.rtg.usage.UsageMetric;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.DiagnosticEvent;
import com.rtg.util.diagnostic.DiagnosticListener;
import com.rtg.util.diagnostic.ListenerType;
import com.rtg.util.test.FileHelper;

import junit.framework.Assert;
import junit.framework.TestCase;

/**
 */
public class NgsTaskFunctionalTest extends TestCase {

  protected File mDir;

  @Override
  public void setUp() throws Exception {
    Diagnostic.setLogStream();
    mDir = FileHelper.createTempDirectory();
  }

  @Override
  public void tearDown() {
    assertTrue(FileHelper.deleteAll(mDir));
    mDir = null;
  }

  static NgsTask getNgs(final NgsParams ngs)  {
    return new NgsTask(ngs, null, new UsageMetric());
  }

  protected NgsParams getParams(final OutputStream out, final NgsMaskParams mask, final ParamsParams restParams, int numThreads) throws IOException {
    return getParams(out, mask, restParams, ListenerType.NULL, OutputFilter.NONE, 10, numThreads);
  }

  NgsParams getParams(OutputStream out, NgsMaskParams mask, ParamsParams restParams, ListenerType listener, OutputFilter filter, int topN, int numThreads) throws IOException {
    final File subjectsDir = FileHelper.createTempDirectory(mDir);
    final File queriesDir = FileHelper.createTempDirectory(mDir);
    final File hitsDir = FileHelper.createTempDirectory(mDir);
    //System.err.println("hitsDir=" + hitsDir);
    ReaderTestUtils.getReaderDNA(restParams.mSubjects, subjectsDir, null).close();
    ReaderTestUtils.getReaderDNA(restParams.mQueries, queriesDir, null).close();
    final SequenceParams subjectParams = SequenceParams.builder().directory(subjectsDir).mode(SequenceMode.UNIDIRECTIONAL).create();
    final SequenceParams queryParams = SequenceParams.builder().directory(queriesDir).create();
    final NgsFilterParams filterParams = NgsFilterParams.builder().outputFilter(filter).zip(restParams.mZip).topN(topN).errorLimit(restParams.mErrorLimit).create();
    final NgsOutputParams outputParams = new OverriddenNgsOutputParams(NgsTestUtils.OverriddenNgsOutputParams.builder().outStream(out).progress(restParams.mProgress).outputDir(new File(hitsDir, "log")).filterParams(filterParams));
    return NgsParams.builder().buildFirstParams(subjectParams).searchParams(queryParams).outputParams(outputParams).maskParams(mask).listeners(Collections.singleton(listener)).compressHashes(false).numberThreads(numThreads).create();
  }

  protected void check(final NgsMaskParams mask, final String subjects, final String queries, final String expected, final Long usageExp) throws Exception {
    check(mask, subjects, queries, expected, null, 1, usageExp);
  }

  void check(final NgsMaskParams mask, final String subjects, final String queries, final String expected, final String[] logExp, final int numThreads, final Long usageExp) throws Exception {
    final ByteArrayOutputStream ba = new ByteArrayOutputStream();
    final PrintStream pr = new PrintStream(ba);
    Diagnostic.setLogStream(pr);
    final ByteArrayOutputStream out = new ByteArrayOutputStream();
    final NgsTask ngs;
    try {
      final NgsParams params = getParams(out, mask, new ParamsParams(subjects, queries, MapFlags.MAX_SCORE, false, false), numThreads);
      ngs = execNgs(params);
    } finally {
      out.close();
      pr.close();
      Diagnostic.setLogStream();
    }
    if (usageExp != null) {
      assertEquals((long) usageExp, ngs.usage());
    }
    final String outString = out.toString();
    if (!expected.equals("")) {
      final String sortede = TestUtils.sortLines(NgsTestUtils.HEADER + expected);
      final String sorteda = TestUtils.sortLines(outString);
      Assert.assertEquals(sortede, sorteda);
    } else {
      Assert.assertEquals("", out.toString());
    }
    if (logExp != null) {
      final String logs = ba.toString();
      //System.err.println(logs);
      TestUtils.containsAll(logs, logExp);
    }
  }


  static NgsTask execNgs(final NgsParams params) throws Exception {
    final NgsTask ngs = getNgs(params);
    final NgsParams par = ngs.parameters();
    assertEquals(params, par);
    ngs.exec();
    params.close();
    return ngs;
  }

  private static final String SEQ_DNA_ODD_S = "" + ">bar" + LS  + "acgt" + LS + ">baaa" + LS + "aaaa";

  private static final String SEQ_DNA_ODD_Q = "" + ">fooo" + LS + "acgt" + LS + ">fiii" + LS + "aaaa";

  public void testOdd() throws Exception {
    //Diagnostic.setLogStream(System.err);
    check(new NgsMaskParamsExplicit("SplitL4w4s0e0"), SEQ_DNA_ODD_S, SEQ_DNA_ODD_Q, "fooo\tF\t0\t1\t0\t0" + LS + "fooo\tR\t0\t1\t0\t0" + LS + "fiii\tF\t1\t1\t0\t0" + LS, null);
  }

  public void testOddg() throws Exception {
    check(new NgsMaskParamsGeneral(4, 0, 0, 1), SEQ_DNA_ODD_S, SEQ_DNA_ODD_Q, "fooo\tF\t0\t1\t0\t0" + LS + "fooo\tR\t0\t1\t0\t0" + LS + "fiii\tF\t1\t1\t0\t0" + LS, null);
  }

  private static final String SEQ_DNA_A = ">x" + LS + "actg" + LS;

  private static final String SEQ_DNA_A2 = ">u" + LS + "actg" + LS + ">v" + LS + "antg" + LS;

  private static final String SEQ_DNA_A3 = ">u" + LS + "actg" + LS + ">v" + LS + "agtg" + LS;

  public void testA1() throws Exception {
    check(new NgsMaskParamsExplicit("SplitL4w4s0e0"), SEQ_DNA_A, SEQ_DNA_A, "x\tF\t0\t1\t0\t0" + LS, null);
  }

  public void testA1g() throws Exception {
    check(new NgsMaskParamsGeneral(4, 0, 0, 1), SEQ_DNA_A, SEQ_DNA_A, "x\tF\t0\t1\t0\t0" + LS, null);
  }

  public void testA2() throws Exception {
    check(new NgsMaskParamsExplicit("SplitL4w2s1e1"), SEQ_DNA_A, SEQ_DNA_A, "x\tF\t0\t1\t0\t0" + LS, null);
  }

  public void testA2g() throws Exception {
    check(new NgsMaskParamsGeneral(2, 1, 1, 1), SEQ_DNA_A, SEQ_DNA_A, "x\tF\t0\t1\t0\t0" + LS, null);
  }

  public void testA3() throws Exception {
    check(new NgsMaskParamsExplicit("SplitL4w2s1e1b"), SEQ_DNA_A2, SEQ_DNA_A, "x\tF\t0\t1\t0\t0" + LS + "x\tF\t1\t1\t1\t1" + LS, null);
  }

  public void testA3g() throws Exception {
    check(new NgsMaskParamsGeneral(2, 1, 1, 1), SEQ_DNA_A2, SEQ_DNA_A, "x\tF\t0\t1\t0\t0" + LS + "x\tF\t1\t1\t1\t1" + LS, null);
  }

  public void testA4() throws Exception {
    check(new NgsMaskParamsExplicit("SplitL4w2s1e1b"), SEQ_DNA_A3, SEQ_DNA_A, "x\tF\t0\t1\t0\t0" + LS + "x\tF\t1\t1\t1\t1" + LS, null);
  }

  public void testA4g() throws Exception {
    check(new NgsMaskParamsGeneral(2, 1, 1, 1), SEQ_DNA_A3, SEQ_DNA_A, "x\tF\t0\t1\t0\t0" + LS + "x\tF\t1\t1\t1\t1" + LS, null);
  }

  public void testA4errorLimit() throws Exception {
    final NgsParams params = getParams(new ByteArrayOutputStream(), new NgsMaskParamsExplicit("SplitL4w2s1e1b"), new NgsTestUtils.ParamsParams(SEQ_DNA_A3, SEQ_DNA_A, 0, false, false), 1);
    execNgs(params);
  }

  public void testA4errorLimitg() throws Exception {
    final NgsParams params = getParams(new ByteArrayOutputStream(), new NgsMaskParamsGeneral(2, 1, 1, 1), new NgsTestUtils.ParamsParams(SEQ_DNA_A3, SEQ_DNA_A, 0, false, false), 1);
    execNgs(params);
  }

  /** Remove the leading "a" to show an example of bug766 */
  private static final String BUG1_TEMPLATE = ">template" + LS + "acccccccctcccccccccccccccccccctccccccc" + LS;

  private static final String BUG1_READ = ">r0" + LS + "cccccccctcaccccccccccccccccccct" + LS;

  public void testBug766() throws Exception {

    check(new NgsMaskParamsGeneral(11, 1, 1, 1), BUG1_READ, BUG1_TEMPLATE, "template\tF\t0\t1\t4\t2" + LS, null);
  }


  private static final String SEQ_DNA_B = ">x" + LS + "acgt" + LS;

  private static final String SEQ_DNA_B2 = ">u" + LS + "acgt" + LS + ">v" + LS + "acgt" + LS;

  public void testB1() throws Exception {
    check(new NgsMaskParamsExplicit("SplitL4w4s0e0"), SEQ_DNA_B, SEQ_DNA_B, "x\tF\t0\t1\t0\t0" + LS + "x\tR\t0\t1\t0\t0" + LS, null);
  }

  public void testB1g() throws Exception {
    check(new NgsMaskParamsGeneral(4, 0, 0, 1), SEQ_DNA_B, SEQ_DNA_B, "x\tF\t0\t1\t0\t0" + LS + "x\tR\t0\t1\t0\t0" + LS, null);
  }

  private static final String READ_WITHN = ">read" + LS + "acgnnnn" + LS;

  private static final String TEMP2_FORWITHN = ""
      + ">temp1" + LS + "acggtac" + LS
      + ">temp2" + LS + "acgaaaa" + LS;

  public void testB2gwithNs() throws Exception {
    check(new NgsMaskParamsGeneral(4, 0, 0, 1), READ_WITHN, TEMP2_FORWITHN, "temp1\tF\t0\t1\t0\t0" + LS, null);
  }

  public void testB2() throws Exception {
    check(new NgsMaskParamsExplicit("SplitL4w2s1e1"), SEQ_DNA_B, SEQ_DNA_B, "x\tF\t0\t1\t0\t0" + LS + "x\tR\t0\t1\t0\t0" + LS, null);
  }

  public void testB2g() throws Exception {
    check(new NgsMaskParamsGeneral(2, 1, 1, 1), SEQ_DNA_B, SEQ_DNA_B, "x\tF\t0\t1\t0\t0" + LS + "x\tR\t0\t1\t0\t0" + LS, null);
  }

  public void testB3() throws Exception {
    check(new NgsMaskParamsExplicit("SplitL4w2s1e1b"), SEQ_DNA_B2, SEQ_DNA_B, "x\tF\t0\t1\t0\t0" + LS + "x\tF\t1\t1\t0\t0" + LS + "x\tR\t0\t1\t0\t0" + LS + "x\tR\t1\t1\t0\t0" + LS, null);
  }

  public void testB3g() throws Exception {
    check(new NgsMaskParamsGeneral(2, 1, 1, 1), SEQ_DNA_B2, SEQ_DNA_B, "x\tF\t0\t1\t0\t0" + LS + "x\tF\t1\t1\t0\t0" + LS + "x\tR\t0\t1\t0\t0" + LS + "x\tR\t1\t1\t0\t0" + LS, null);
  }

  //These carefully test indel handling including variable length indels
  //genome
  private static final String SEQ_DNA_C_GE = ">x" + LS + "aaccggttaacc" + "ggaatt" + LS;

  //reads with carefully selected substitutions and indels
  //The capitalized nt are the ones that will match
  private static final String SEQ_DNA_C_RE = ""
      + ">r0" + LS + "gcTTatCCacTT" + LS
      + ">r1" + LS + "gcTTatGGacTT" + LS
      + ">r2" + LS + "gcTTatCGacTT" + LS
      ;
  private void checkC(final int indel, final int l, final String[] exp, final String[] nex) throws Exception {
    final ByteArrayOutputStream out = new ByteArrayOutputStream();
    try {
      final NgsParams params = getParams(out, new NgsMaskParamsGeneral(6, 3, indel, l), new NgsTestUtils.ParamsParams(SEQ_DNA_C_RE, SEQ_DNA_C_GE, MapFlags.MAX_SCORE, false, false), 1);
      execNgs(params);
    } finally {
      out.close();
    }
    final String actual = out.toString();
    //System.err.println(actual);
    for (int i = 0; i < exp.length; i++) {
      final String str = exp[i];
      if (!actual.contains(str)) {
        fail("String not found:" + i + LS + str + LS + actual);
      }
    }
    for (int i = 0; i < nex.length; i++) {
      final String str = nex[i];
      if (actual.contains(str)) {
        fail("String found:" + i + LS + str + LS + actual);
      }
    }
  }


  public void testC1() throws Exception {
    final String[] exp = {
    };
    final String[] nex = {
    };
    checkC(0, 1, exp, nex);
  }

  public void testC2() throws Exception {
    final String[] exp = {
        "x\tF\t0\t7",
        "x\tF\t1\t7",
        "x\tF\t2\t7",
    };
    final String[] nex = {
    };
    checkC(2, 1, exp, nex);
  }

  public void testC3() throws Exception {
    final String[] exp = {
        "x\tF\t0\t7",
        "x\tF\t1\t7",
    };
    final String[] nex = {
        "x\tF\t2\t7",
    };
    checkC(1, 2, exp, nex);
  }

  private static final String SEQ_DNA_U1 = ">x" + LS + "actg" + LS;

  private static final String SEQ_DNA_U2 = ">u" + LS + "actg" + LS;

  public void testU() throws Exception {
    check(new NgsMaskParamsExplicit("SplitL4w2s1e1"), SEQ_DNA_U2, SEQ_DNA_U1, "x\tF\t0\t1\t0\t0" + LS, null);
  }

  public void testUg() throws Exception {
    check(new NgsMaskParamsGeneral(2, 1, 1, 1), SEQ_DNA_U2, SEQ_DNA_U1, "x\tF\t0\t1\t0\t0" + LS, null);
  }

  private static final String SEQ_DNA_V1 = ">x" + LS + "ttactg" + LS;

  private static final String SEQ_DNA_V2 = ">u" + LS + "actg" + LS;

  /** tzero != tstart case */
  public void testV() throws Exception {
    check(new NgsMaskParamsExplicit("SplitL4w4s0e0"), SEQ_DNA_V2, SEQ_DNA_V1, "x\tF\t0\t3\t0\t0" + LS, null);
  }

  public void testVg() throws Exception {
    check(new NgsMaskParamsGeneral(4, 0, 0, 1), SEQ_DNA_V2, SEQ_DNA_V1, "x\tF\t0\t3\t0\t0" + LS, null);
  }

  private static final String SEQ_DNA_W1 = ">x" + LS + "actgactg" + LS;

  private static final String SEQ_DNA_W2 = ">u" + LS + "actg" + LS;

  /** tzero != 1 case */
  public void testW() throws Exception {
    check(new NgsMaskParamsExplicit("SplitL4w4s0e0"), SEQ_DNA_W2, SEQ_DNA_W1, "x\tF\t0\t1\t0\t0" + LS + "x\tF\t0\t5\t0\t0" + LS, null);
  }

  /** tzero != 1 case */
  public void testWg() throws Exception {
    check(new NgsMaskParamsGeneral(4, 0, 0, 1), SEQ_DNA_W2, SEQ_DNA_W1, "x\tF\t0\t1\t0\t0" + LS + "x\tF\t0\t5\t0\t0" + LS, null);
  }

  /** Subject sequence.  */
  public static final String SEQ_X1 = ">x0" + LS + "actgactg" + LS + ">x1" + LS + "actgactg" + LS + ">x2" + LS + "actgactgactg" + LS + ">x3" + LS + "actg" + LS + ">x4" + LS + "actgactg" + LS + ">x5" + LS + "actgactg" + LS + ">x6" + LS + "actg" + LS + ">x7" + LS  + "actgactgactg" + LS + ">x8" + LS + "actgactgactg" + LS + ">x9" + LS + "actgactgactgactgactgactgactgactgactg" + LS;
  /** Query sequence.  */
  static final String SEQ_X2 = ">u" + LS + "actg" + LS;

  public void testX() throws Exception {
    final DiagnosticListener l = new DiagnosticListener() {
      @Override
      public void handleDiagnosticEvent(final DiagnosticEvent<?> event) {
        if (event.getMessage().equals("Preread")) {
          return;
        }
      }

      @Override
      public void close() {
      }
    };
    Diagnostic.addListener(l);
    try {
      check(new NgsMaskParamsExplicit("SplitL36w18s2e1"), SEQ_X1, SEQ_X2, "", null);
    } finally {
      Diagnostic.removeListener(l);
    }
  }

  /** Subject sequence.  */
  static final String BUG_READS = ">x0" + LS + "actg" + LS + ">x1" + LS + "actg" + LS + ">x2" + LS + "gctg" + LS + ">x3" + LS + "gctg" + LS + ">x4" + LS + "gctg" + LS + ">x5" + LS + "cctg" + LS + ">x6" + LS + "cctg" + LS + ">x7" + LS + "cctg" + LS + ">x8" + LS + "cctg" + LS + ">x9" + LS + "acgg" + LS;

  /** Query sequence.  */
  static final String BUG_TEMPLATE = ">bug" + LS + "actg" + LS;

  static final String EXPECTED_BUG = "bug\tF\t0\t1\t0\t0" + LS + "bug\tF\t1\t1\t0\t0" + LS + "bug\tF\t9\t1\t1\t1" + LS;

  /** Picks up a bug seen in the wild. */
  public void testBug() throws Exception {
    checkBug(new NgsMaskParamsExplicit("SplitL4w2s1e1b"));
  }

  /** Picks up a bug seen in the wild. */
  public void testBugg() throws Exception {
    checkBug(new NgsMaskParamsGeneral(2, 1, 1, 1));
  }

  private void checkBug(final NgsMaskParams mask) throws Exception {
    final File subjectsDir = FileHelper.createTempDirectory(mDir);
    final File queriesDir = FileHelper.createTempDirectory(mDir);
    final File hitsDir = FileHelper.createTempDirectory(mDir);
    //System.err.println("hitsDir=" + hitsDir);
    final ByteArrayOutputStream out = new ByteArrayOutputStream();
    final Appendable err = new StringWriter();
    final SequencesReader sr = ReaderTestUtils.getReaderDNA(BUG_READS, subjectsDir, null);
    try {
      final SequencesReader qr = ReaderTestUtils.getReaderDNA(BUG_TEMPLATE, queriesDir, null);
      try {
        final SequenceParams subjectParams = SequenceParams.builder().directory(subjectsDir).mode(SequenceMode.UNIDIRECTIONAL).create();
        final SequenceParams queryParams = SequenceParams.builder().directory(queriesDir).create();
        final NgsFilterParams filterParams = NgsFilterParams.builder().outputFilter(OutputFilter.NONE).topN(10).errorLimit(5).create();
        final NgsOutputParams outputParams = new NgsTestUtils.OverriddenNgsOutputParams(NgsTestUtils.OverriddenNgsOutputParams.builder().outStream(out).progress(false).outputDir(new File(hitsDir, "log")).filterParams(filterParams));
        final NgsParams params = NgsParams.builder().hashCountThreshold(5).buildFirstParams(subjectParams).searchParams(queryParams).outputParams(outputParams).maskParams(mask).create(); //less than the number of reads
        final NgsTask ngs = getNgs(params);
        final NgsParams par = ngs.parameters();
        assertEquals(params, par);
        ngs.exec();

        assertEquals(NgsTestUtils.HEADER + EXPECTED_BUG, out.toString() + err.toString());
        queryParams.close();
        subjectParams.close();
      } finally {
        qr.close();
      }
    } finally {
      sr.close();
    }
  }

  private static final String BUG2_BUILD = ">x" + LS + "TGGTCAGGTTCTGCGTGACGAGCTACTGCGTCGCGGCTGGATCCACCAATT" + LS;

  private static final String BUG2_SEARCH = ">u" + LS + "AGCATCACATCTCTGGTCAGGTTCTGCGTGACGAGCTACTGCGTCGCGGCTGGATCCACCAATTACGTCAACGCAGCGACGGTC" + LS
      ;

  private static final String EXPECTED_BUG2_A = "x\tF\t0\t1\t0\t0" + LS;

  private static final String EXPECTED_BUG2_B = "u\tF\t0\t14\t0\t0" + LS;

  public void testBug2a() throws Exception {
    check(new NgsMaskParamsExplicit("MaskL51w12s3e3"), BUG2_BUILD, BUG2_BUILD, EXPECTED_BUG2_A, null);
  }

  public void testBug2ag() throws Exception {
    check(new NgsMaskParamsGeneral(12, 3, 3, 1), BUG2_BUILD, BUG2_BUILD, EXPECTED_BUG2_A, null);
  }

  public void testBug2b() throws Exception {
    check(new NgsMaskParamsExplicit("MaskL51w12s3e3"), BUG2_BUILD, BUG2_SEARCH, EXPECTED_BUG2_B, null);
  }

  public void testBug2bg() throws Exception {
    check(new NgsMaskParamsGeneral(12, 3, 3, 1), BUG2_BUILD, BUG2_SEARCH, EXPECTED_BUG2_B, null);
  }

  protected static final String SEQ_DNA_Y1 = ">x" + LS + "actg" + LS;

  protected static final String SEQ_DNA_Y2 = ">u" + LS + "actg" + LS + ">v" + LS + "actg" + LS;

  protected static final String EXPECTED_Y_Z = "x\tF\t0\t1\t0\t0" + LS + "x\tF\t1\t1\t0\t0" + LS;

  /**
   * Should be the same as Z
   */
  public void testY() throws Exception {
    check(new NgsMaskParamsExplicit("SplitL4w4s0e0"), SEQ_DNA_Y2, SEQ_DNA_Y1, EXPECTED_Y_Z, null);
  }

  /**
   * Should be the same as Z
   */
  public void testYg() throws Exception {
    check(new NgsMaskParamsGeneral(4, 0, 0, 1), SEQ_DNA_Y2, SEQ_DNA_Y1, EXPECTED_Y_Z, null);
  }

  private static final String READS_36 = ">r36" + LS  + "aaaaaaaaaccccccccctttttttttggggggggg" + LS;

  private static final String TEMPLATE_36_0 = ">t36" + LS + "tttttttttttttttttttttttttttggggggggg" + LS;

  private static final String TEMPLATE_36_1 = ">t36" + LS + "cccccccccccccccccctttttttttttttttttt" + LS;

  private static final String TEMPLATE_36_2 = ">t36" + LS + "aaaaaaaaaccccccccccccccccccccccccccc" + LS;

  private static final String TEMPLATE_36_3 = ">t36" + LS + "cccccccccccccccccccccccccccggggggggg" + LS;

  private static final String TEMPLATE_36_4 = ">t36" + LS + "aaaaaaaaaaaaaaaaaatttttttttttttttttt" + LS;

  private static final String TEMPLATE_36_5 = ">t36" + LS + "aaaaaaaaaaaaaaaaaagggggggggggggggggg" + LS;

  public void test36a() throws Exception {
    check(new NgsMaskParamsExplicit("SplitL36w18s2e1"), READS_36, TEMPLATE_36_0, "t36\tF\t0\t1\t18\t18" + LS, null);
  }

  public void test36ag() throws Exception {
    check(new NgsMaskParamsGeneral(18, 2, 1, 1), READS_36, TEMPLATE_36_0, "t36\tF\t0\t1\t18\t18" + LS, null);
  }

  public void test36b() throws Exception {
    check(new NgsMaskParamsExplicit("SplitL36w18s2e1"), READS_36, TEMPLATE_36_1, "" + "t36\tF\t0\t1\t18\t18" + LS + "t36\tR\t0\t1\t18\t18" + LS, null);
  }

  public void test36bg() throws Exception {
    check(new NgsMaskParamsGeneral(18, 2, 1, 1), READS_36, TEMPLATE_36_1, "" + "t36\tF\t0\t1\t18\t18" + LS + "t36\tR\t0\t1\t18\t18" + LS, null);
  }

  public void test36c() throws Exception {
    check(new NgsMaskParamsExplicit("SplitL36w18s2e1"), READS_36, TEMPLATE_36_2, "t36\tF\t0\t1\t18\t18" + LS, null);
  }

  public void test36cg() throws Exception {
    check(new NgsMaskParamsGeneral(18, 2, 1, 1), READS_36, TEMPLATE_36_2, "t36\tF\t0\t1\t18\t18" + LS, null);
  }

  public void test36d() throws Exception {
    check(new NgsMaskParamsExplicit("SplitL36w18s2e1"), READS_36, TEMPLATE_36_3, "t36\tF\t0\t1\t18\t18" + LS, null);
  }

  public void test36dg() throws Exception {
    check(new NgsMaskParamsGeneral(18, 2, 1, 1), READS_36, TEMPLATE_36_3, "t36\tF\t0\t1\t18\t18" + LS, null);
  }

  public void test36e() throws Exception {
    check(new NgsMaskParamsExplicit("SplitL36w18s2e1"), READS_36, TEMPLATE_36_4, "" + "t36\tF\t0\t1\t18\t17" + LS + "t36\tR\t0\t1\t18\t17" + LS, null);
  }

  public void test36eg() throws Exception {
    check(new NgsMaskParamsGeneral(18, 2, 1, 1), READS_36, TEMPLATE_36_4, "" + "t36\tF\t0\t1\t18\t17" + LS + "t36\tR\t0\t1\t18\t17" + LS, null);
  }

  public void test36f() throws Exception {
    check(new NgsMaskParamsExplicit("SplitL36w18s2e1"), READS_36, TEMPLATE_36_5, "" + "t36\tF\t0\t1\t18\t18" + LS + "t36\tR\t0\t1\t18\t18" + LS, null);
  }

  public void test36fg() throws Exception {
    check(new NgsMaskParamsGeneral(18, 2, 1, 1), READS_36, TEMPLATE_36_5, "" + "t36\tF\t0\t1\t18\t18" + LS + "t36\tR\t0\t1\t18\t18" + LS, null);
  }

  private static final String READ_32_DEL_SEQ = "c" + "ttgcagtatagcg tgctagcat catgagcg" + "a";
  //0     1234567890123 456789012 34567890     1
  private static final String READ_32_DEL = ""
      + ">read" + LS
      + READ_32_DEL_SEQ + LS;
  private static final String READ_32_DEL_REVERSE = ""
      + ">read" + LS
      + DnaUtils.reverseComplement(READ_32_DEL_SEQ) + LS;

  private static final String TEMPLATE_32_DEL_SEQ = "atc ttgcagtatagcg a tgctagcat t catgagcg atcgatctatcgattattgca";
  private static final String TEMPLATE_32_DEL = ""
      + ">template" + LS
      + TEMPLATE_32_DEL_SEQ + LS;

  public void test32Del12() throws Exception {
    //each deletion seems to give +1 to start position, so the actual start of 3 becomes 5
    check(new NgsMaskParamsGeneral(12, 2, 2, 1), READ_32_DEL, TEMPLATE_32_DEL, "template\tF\t0\t5\t19\t11" + LS, null);
    check(new NgsMaskParamsGeneral(12, 2, 2, 1), READ_32_DEL_REVERSE, TEMPLATE_32_DEL,
        "template\tR\t0\t3\t16\t8" + LS, null);
  }

  public void test32Del16() throws Exception {
    //each deletion seems to give +1 to start position, so the actual start of 3 becomes 5
    check(new NgsMaskParamsGeneral(16, 2, 2, 1), READ_32_DEL, TEMPLATE_32_DEL, "template\tF\t0\t5\t19\t11" + LS, null);
    check(new NgsMaskParamsGeneral(16, 2, 2, 1), READ_32_DEL_REVERSE, TEMPLATE_32_DEL,
        "template\tR\t0\t3\t16\t8" + LS, null);
  }

  private static final String BLOCK_1 = "acgtccggt";
  private static final String BLOCS_1 = "acctccggt";
  private static final String BLOCK_2 = "acggtagcg";
  private static final String BLOCK_3 = "acgtccggt";
  private static final String BLOCK_4 = "gatacagat";

  private static final String READS_BLOCK = ""
      + ">t36" + LS
      + BLOCK_1 + BLOCK_2 + BLOCK_3 + BLOCK_4 + LS;

  private static final String JUNK_1 =  "ggaaccccc";
  private static final String JUNK_2 =  "gggggcccc";
  private static final String PRE_TEMPLATE = ">t36" + LS + "ctgtac";

  public void test36ba() throws Exception {
    final String tem = PRE_TEMPLATE + JUNK_1 + JUNK_2 + BLOCK_3 + BLOCK_4 + LS;
    check(new NgsMaskParamsExplicit("SplitL36w18s2e1"), READS_BLOCK, tem, "t36\tF\t0\t7\t13\t11" + LS, null);
  }

  public void test36bag() throws Exception {
    final String tem = PRE_TEMPLATE + JUNK_1 + JUNK_2 + BLOCK_3 + BLOCK_4 + LS;
    check(new NgsMaskParamsGeneral(18, 2, 1, 1), READS_BLOCK, tem, "t36\tF\t0\t7\t13\t11" + LS, null);
  }

  public void test36bb() throws Exception {
    final String tem = PRE_TEMPLATE + JUNK_1 + BLOCK_2 + BLOCK_3 + JUNK_2 + LS;
    check(new NgsMaskParamsExplicit("SplitL36w18s2e1"), READS_BLOCK, tem, "t36\tF\t0\t7\t15\t13" + LS, null);
  }

  public void test36bbg() throws Exception {
    final String tem = PRE_TEMPLATE + JUNK_1 + BLOCK_2 + BLOCK_3 + JUNK_2 + LS;
    check(new NgsMaskParamsGeneral(18, 2, 1, 1), READS_BLOCK, tem, "t36\tF\t0\t7\t15\t13" + LS, null);
  }

  public void test36bc() throws Exception {
    final String tem = PRE_TEMPLATE + BLOCK_1 + BLOCK_2 + JUNK_1 + JUNK_2 + LS;
    check(new NgsMaskParamsExplicit("SplitL36w18s2e1"), READS_BLOCK, tem, "t36\tF\t0\t7\t15\t12" + LS, null);
  }

  public void test36bcg() throws Exception {
    final String tem = PRE_TEMPLATE + BLOCK_1 + BLOCK_2 + JUNK_1 + JUNK_2 + LS;
    check(new NgsMaskParamsGeneral(18, 2, 1, 1), READS_BLOCK, tem, "t36\tF\t0\t7\t15\t12" + LS, null);
  }

  public void test36bd() throws Exception {
    final String tem = PRE_TEMPLATE + JUNK_1 + BLOCK_2 + JUNK_2 + BLOCK_4 + LS;
    check(new NgsMaskParamsExplicit("SplitL36w18s2e1"), READS_BLOCK, tem, "t36\tF\t0\t7\t14\t11" + LS, null);
  }

  public void test36bdg() throws Exception {
    final String tem = PRE_TEMPLATE + JUNK_1 + BLOCK_2 + JUNK_2 + BLOCK_4 + LS;
    check(new NgsMaskParamsGeneral(18, 2, 1, 1), READS_BLOCK, tem, "t36\tF\t0\t7\t14\t11" + LS, null);
  }

  public void test36be() throws Exception {
    final String tem = PRE_TEMPLATE + BLOCK_1 + JUNK_2 + BLOCK_3 + JUNK_1 + LS;
    check(new NgsMaskParamsExplicit("SplitL36w18s2e1"), READS_BLOCK, tem, "t36\tF\t0\t7\t12\t8" + LS, null);
  }

  public void test36beg() throws Exception {
    final String tem = PRE_TEMPLATE + BLOCK_1 + JUNK_2 + BLOCK_3 + JUNK_1 + LS;
    check(new NgsMaskParamsGeneral(18, 2, 1, 1), READS_BLOCK, tem, "t36\tF\t0\t7\t12\t8" + LS, null);
  }

  public void test36bf() throws Exception {
    final String tem = PRE_TEMPLATE + BLOCK_1 + JUNK_2 + JUNK_1 + BLOCK_4 + LS;
    check(new NgsMaskParamsExplicit("SplitL36w18s2e1"), READS_BLOCK, tem, "t36\tF\t0\t7\t13\t10" + LS, null);
  }

  public void test36bfg() throws Exception {
    final String tem = PRE_TEMPLATE + BLOCK_1 + JUNK_2 + JUNK_1 + BLOCK_4 + LS;
    check(new NgsMaskParamsGeneral(18, 2, 1, 1), READS_BLOCK, tem, "t36\tF\t0\t7\t13\t10" + LS, null);
  }

  private static String reverseComplement(final String nt) {
    final StringBuilder sb = new StringBuilder();
    for (int j = nt.length() - 1; j >= 0; j--) {
      final char c = nt.charAt(j);
      final char d;
      switch (c) {
        case 'a':
        case 'A': d = 't'; break;
        case 'c':
        case 'C':
          d = 'g'; break;
        case 'g':
        case 'G':
          d = 'c'; break;
        case 't':
        case 'T':
          d = 'a'; break;
        default: throw new RuntimeException();
      }
      sb.append(d);
    }
    return sb.toString();
  }

  private static final String POST = "cac";

  private static final String READ_RC = ""
      + ">r1" + LS
      + "cgat" + LS;

  private static final String TEMP_RC = ""
      + ">t1" + LS
      + "atcg" + LS;

  public void testRC() throws Exception {
    check(new NgsMaskParamsExplicit("SplitL4w4s0e0"), READ_RC, TEMP_RC, "t1\tR\t0\t1\t0\t0" + LS, 4L);
  }

  public void testRCg() throws Exception {
    check(new NgsMaskParamsGeneral(4, 0, 0, 1), READ_RC, TEMP_RC, "t1\tR\t0\t1\t0\t0" + LS, null);
  }

  private static final String READ_BUG = ">r1" + LS + "TTACGACCAATGAGCGTCGCCAGCCCCATCCCCAGC" + LS;

  private static final String TEMP_BUG = ">t1" + LS + "CGCTGGGGATGGGGCTGGCGACGCTGATTGGTCGTAAAA" + LS;

  public void testRCBug() throws Exception {
    check(new NgsMaskParamsExplicit("SplitL36w18s2e1"), READ_BUG, TEMP_BUG, "t1\tR\t0\t2\t1\t1" + LS, 36L);
  }

  public void testRCBugg() throws Exception {
    check(new NgsMaskParamsExplicit("SplitL36w18s2e1"), READ_BUG, TEMP_BUG, "t1\tR\t0\t2\t1\t1" + LS, null);
  }

  public void test36rc() throws Exception {
    final String seq = BLOCK_1 + BLOCK_2 + BLOCK_3 + BLOCK_4;
    final String tem = PRE_TEMPLATE + reverseComplement(seq) + POST + LS;
    check(new NgsMaskParamsExplicit("SplitL36w18s2e1"), READS_BLOCK, tem, "t36\tR\t0\t7\t0\t0" + LS, (long) seq.length());
  }

  public void test36rcg() throws Exception {
    final String tem = PRE_TEMPLATE + reverseComplement(BLOCK_1 + BLOCK_2 + BLOCK_3 + BLOCK_4) + POST + LS;
    check(new NgsMaskParamsGeneral(18, 2, 1, 1), READS_BLOCK, tem, "t36\tR\t0\t7\t0\t0" + LS, null);
  }

  public void test36rcs() throws Exception {
    final String tem = PRE_TEMPLATE + reverseComplement(BLOCS_1 + BLOCK_2 + BLOCK_3 + BLOCK_4) + POST + LS;
    check(new NgsMaskParamsExplicit("SplitL36w18s2e1"), READS_BLOCK, tem, "t36\tR\t0\t7\t1\t1" + LS, null);
  }

  public void test36rcsg() throws Exception {
    final String tem = PRE_TEMPLATE + reverseComplement(BLOCS_1 + BLOCK_2 + BLOCK_3 + BLOCK_4) + POST + LS;
    check(new NgsMaskParamsGeneral(18, 2, 1, 1), READS_BLOCK, tem, "t36\tR\t0\t7\t1\t1" + LS, null);
  }

  public void test36rca() throws Exception {
    final String seq = JUNK_1 + JUNK_2 + BLOCK_3 + BLOCK_4;
    final String tem = PRE_TEMPLATE + reverseComplement(seq) + POST + LS;
    //TODO check this more carefully
    check(new NgsMaskParamsExplicit("SplitL36w18s2e1"), READS_BLOCK, tem, "t36\tR\t0\t7\t13\t11" + LS, (long) seq.length());
  }

  public void test36rcag() throws Exception {
    final String tem = PRE_TEMPLATE + reverseComplement(JUNK_1 + JUNK_2 + BLOCK_3 + BLOCK_4) + POST + LS;
    check(new NgsMaskParamsGeneral(18, 2, 1, 1), READS_BLOCK, tem, "t36\tR\t0\t7\t13\t11" + LS, null);
  }

  public void test36rcb() throws Exception {
    final String tem = PRE_TEMPLATE + reverseComplement(JUNK_1 + BLOCK_2 + BLOCK_3 + JUNK_2) + POST + LS;
    //TODO check this more carefully
    check(new NgsMaskParamsExplicit("SplitL36w18s2e1"), READS_BLOCK, tem, "t36\tR\t0\t7\t15\t13" + LS, null);
  }

  public void test36rcbg() throws Exception {
    final String tem = PRE_TEMPLATE + reverseComplement(JUNK_1 + BLOCK_2 + BLOCK_3 + JUNK_2) + POST + LS;
    check(new NgsMaskParamsGeneral(18, 2, 1, 1), READS_BLOCK, tem, "t36\tR\t0\t7\t15\t13" + LS, null);
  }

  public void test36rcc() throws Exception {
    final String tem = PRE_TEMPLATE + reverseComplement(BLOCK_1 + BLOCK_2 + JUNK_1 + JUNK_2) + POST + LS;
    check(new NgsMaskParamsExplicit("SplitL36w18s2e1"), READS_BLOCK, tem, "t36\tR\t0\t7\t15\t12" + LS, null);
  }

  public void test36rccg() throws Exception {
    final String tem = PRE_TEMPLATE + reverseComplement(BLOCK_1 + BLOCK_2 + JUNK_1 + JUNK_2) + POST + LS;
    check(new NgsMaskParamsGeneral(18, 2, 1, 1), READS_BLOCK, tem, "t36\tR\t0\t7\t15\t12" + LS, null);
  }

  public void test36rcd() throws Exception {
    final String tem = PRE_TEMPLATE + reverseComplement(JUNK_1 + BLOCK_2 + JUNK_2 + BLOCK_4) + POST + LS;
    //TODO check this more carefully
    check(new NgsMaskParamsExplicit("SplitL36w18s2e1"), READS_BLOCK, tem, "t36\tR\t0\t7\t14\t11" + LS, 36L);
  }

  public void test36rcdg() throws Exception {
    final String tem = PRE_TEMPLATE + reverseComplement(JUNK_1 + BLOCK_2 + JUNK_2 + BLOCK_4) + POST + LS;
    check(new NgsMaskParamsGeneral(18, 2, 1, 1), READS_BLOCK, tem, "t36\tR\t0\t7\t14\t11" + LS, null);
  }

  public void test36rce() throws Exception {
    final String tem = PRE_TEMPLATE + reverseComplement(BLOCK_1 + JUNK_2 + BLOCK_3 + JUNK_1) + POST + LS;
    check(new NgsMaskParamsExplicit("SplitL36w18s2e1"), READS_BLOCK, tem, "t36\tR\t0\t7\t12\t8" + LS, null);
  }

  public void test36rceg() throws Exception {
    final String tem = PRE_TEMPLATE + reverseComplement(BLOCK_1 + JUNK_2 + BLOCK_3 + JUNK_1) + POST + LS;
    check(new NgsMaskParamsGeneral(18, 2, 1, 1), READS_BLOCK, tem, "t36\tR\t0\t7\t12\t8" + LS, null);
  }

  public void test36rcf() throws Exception {
    final String tem = PRE_TEMPLATE + reverseComplement(BLOCK_1 + JUNK_2 + JUNK_1 + BLOCK_4) + POST + LS;
    check(new NgsMaskParamsExplicit("SplitL36w18s2e1"), READS_BLOCK, tem, "t36\tR\t0\t7\t13\t10" + LS, null);
  }

  public void test36rcfg() throws Exception {
    final String tem = PRE_TEMPLATE + reverseComplement(BLOCK_1 + JUNK_2 + JUNK_1 + BLOCK_4) + POST + LS;
    check(new NgsMaskParamsGeneral(18, 2, 1, 1), READS_BLOCK, tem, "t36\tR\t0\t7\t13\t10" + LS, null);
  }

  static final String READS_BUG = ">read53:gi|48994873|gb|U00096.2|:4415055:S0:I1:D0" + LS + "ACCGTGCGTAATTTTTTATCACGGCTTaTACTTCAT" + LS;

  private static final String TEM_BUG = ">tbug" + LS + "cgttACCGTGCGTAATTTTTTATCACGGCTTTACTTCATcc" + LS;

  /** Bug found in the wild. */
  public void testBug36() throws Exception {
    check(new NgsMaskParamsExplicit("SplitL36w18s2e1"), READS_BUG, TEM_BUG, "tbug\tF\t0\t5\t8\t2" + LS, 36L);
  }

  /** Bug found in the wild. */
  public void testBug36g() throws Exception {
    check(new NgsMaskParamsGeneral(18, 2, 1, 1), READS_BUG, TEM_BUG, "tbug\tF\t0\t5\t8\t2" + LS, null);
  }

  public void testMultiCoreSplitting50() throws Exception {
    final String template = makeReads(">seq@" + LS + "ACGT" + LS, 50);
    final String reads = ""
        + ">read" + LS
        + "ACTG" + LS
        ;
    final String x = ""
        + "seq@" + TAB + "F" + TAB + "0" + TAB + "1" + TAB + "2" + TAB + "2" + LS
        + "seq@" + TAB + "R" + TAB + "0" + TAB + "1" + TAB + "2" + TAB + "2" + LS
        ;
    final String exp = makeReads(x, 50);

    final String[] logExp = {
        "threads=1",
    };
    check(new NgsMaskParamsGeneral(2, 1, 1, 1), reads, template, exp, logExp, 1, 4L);
  }

  static String makeReads(final String str, final int repeat) {
    final StringBuilder sb = new StringBuilder();
    for (int i = 0; i < repeat; i++) {
      sb.append(str.replace("@", "" + i));
    }
    return sb.toString();
  }

  public void testMultiCoreSplitting40() throws Exception {
    final int max = 4 * NgsHashLoopImpl.MAX_SEQUENCES;
    final String template = makeReads(">seq@" + LS + "ACGT" + LS, max);
    final String reads = ""
        + ">read" + LS
        + "ACTG" + LS
        ;
    final String x = ""
        + "seq@" + TAB + "F" + TAB + "0" + TAB + "1" + TAB + "2" + TAB + "2" + LS
        + "seq@" + TAB + "R" + TAB + "0" + TAB + "1" + TAB + "2" + TAB + "2" + LS
        ;
    final String exp = makeReads(x, max);

    final String[] logExp = {
        "threads=1",
    };
    check(new NgsMaskParamsGeneral(2, 1, 1, 1), reads, template, exp, logExp, 1, 4L);
  }

  private static final String RB3 = "GCGACGCTGCCGCTGGTGTAGGGTTCCGCTTTGTTT";

  private static final String READS_BUGE3 = ">read1418:gi|48994873|gb|U00096.2|:3477973R:S3:I0:D0" + LS + RB3 + LS;

  private static final String TB3 = "tAAACAAtGCGGAACCCaACACCAGCaGCAGCGTCGCGg";

  private static final String TEM_BUGE3 = ">tbug" + LS + TB3 + LS;

  /** Bug in e=3 mask. Seen in the wild. */
  public void testBug36e3() throws Exception {
    check(new NgsMaskParamsExplicit("SplitL36w12s3e3"), READS_BUGE3, TEM_BUGE3, "tbug\tR\t0\t2\t3\t3" + LS, (long) RB3.length());
  }

  /** Bug in e=3 mask. Seen in the wild. */
  public void testBug36e3g() throws Exception {
    check(new NgsMaskParamsGeneral(12, 3, 3, 1), READS_BUGE3, TEM_BUGE3, "tbug\tR\t0\t2\t3\t3" + LS,  (long) RB3.length());
  }

  public void testC1Log() throws Exception {
    final ByteArrayOutputStream ba = new ByteArrayOutputStream();
    final PrintStream pr = new PrintStream(ba);
    try {
      Diagnostic.setLogStream(pr);
      try {
        final ByteArrayOutputStream out = new ByteArrayOutputStream();
        final NgsParams params = getParams(out, new NgsMaskParamsExplicit("SplitL4w4s0e0"), new NgsTestUtils.ParamsParams(SEQ_DNA_Y2, SEQ_DNA_Y1, 10/*errorlimit*/, true/*progress*/, false/*zip*/), ListenerType.NULL,  OutputFilter.NONE, 1/*topN*/, 1);
        final NgsTask ngs;
        try {
          ngs = execNgs(params);
        } finally {
          out.close();
        }
        assertEquals(8, ngs.usage());
        assertEquals(NgsTestUtils.HEADER + EXPECTED_Y_Z, out.toString());
        //final File log = params.output().file("log");
        pr.flush();
        final String logs = ba.toString();
        //System.err.println(logs);
        //read counts
        assertTrue(logs, logs.contains(" filter=NONE"));
        //other timers
        //        assertTrue(logs.contains(" Timer Read_delay_subject "));
        if (params.numberThreads() > 1) {
          TestUtils.containsAll(logs, "Thread Search ", " Thread Search parent Terminating ", " Thread Search parent Finished ",
            " Thread Search 0 Scheduling", " Thread Search 0 Start ", " Thread Search 0 Read_done ", " Thread Search 0 Finish ",
            " Thread Search 0R Scheduling", " Thread Search 0R Start ", " Thread Search 0R Finish ");
        }
        //deprecated timers to be removed
        TestUtils.containsAll(logs, " Timer Index_initialization ", " Timer Index_sort ", " Timer Index_pointer ", " Timer Index_position ", " Timer Index_bitVector ");
        //index statistics
        TestUtils.containsAll(logs, "Index[", "] memory performance", "] search performance");
      } finally {
        Diagnostic.setLogStream();
      }
    } finally {
      pr.close();
      ba.close();
    }
  }
}

