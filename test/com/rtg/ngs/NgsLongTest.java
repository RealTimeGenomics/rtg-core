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

import com.rtg.ngs.NgsTestUtils.TestParams;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.LogRecord;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * Tests for Long mapping
 */
public class NgsLongTest extends TestCase {


  public static void main(final String[] args) {
    junit.textui.TestRunner.run(suite());
  }

  public static Test suite() {
    final TestSuite suite = new TestSuite();

    suite.addTestSuite(NgsLongTest.class);
    return suite;
  }

  @Override
  public void setUp() {
    Diagnostic.setLogStream();
  }

  //  private static final String READ_1 = ">read_l/1" + LS
  //  + "TGCCAC" + LS
  //  + ">read_2/2" + LS + "TGCTGT" + LS;
  //12345678901234567890123456789012345678901234567890123456789012

  //  private static final String TEMP_1 = ">template" + LS
  //  + "CAGGCAACTGCCACCTT" + LS
  //  + "GGTTTTTTGCCCCCCTG" + LS;

  //  private static final String EXP_LONG1TOPN = ""
  //    + "template\tF\t0\t9\t0\t0" + LS
  //    + "template\tR\t0\t1\t0\t1" + LS
  //    + "template\tF\t0\t25\t0\t1" + LS
  //    + "#read 0 had 10 results with score-indel 2" + LS
  //    + "#read 1 had 10 results with score-indel 2" + LS;

  //  public void testLong1Topn() {
  //    final int r = 6, w = 2, a = 0, b = 0, c = 1;
  //    final boolean cgl = false;
  //    final int stepSize = 2;
  //    final boolean useLongReads = true;
  //    checkLong(new NgsMaskParamsGeneral(w, a, b, c, cgl),
  //        new NgsTestUtils.TestParams(useLongReads,
  //            stepSize, READ_1, TEMP_1, EXP_LONG1TOPN),
  //            new NgsTestUtils.NgsFilterPartlyParams(OutputFilter.TOPN, 5), false);
  //  }

  //  private static final String EXP_SHORT1TOPN = ""
  //    + "template\tF\t0\t9\t0\t0" + LS
  //    + "template\tF\t0\t3\t2\t2" + LS
  //    + "template\tR\t0\t19\t4\t4" + LS
  //    + "template\tR\t1\t7\t2\t2" + LS
  //    + "template\tF\t1\t15\t4\t4" + LS
  //    + "template\tR\t1\t13\t4\t4" + LS;

  //  public void testShort1Topn() {
  //    final int r = 6, w = 2, a = 0, b = 0, c = 1;
  //    final boolean cgl = false;
  //    NgsTestUtils.checkShort(new NgsMaskParamsGeneral(w, a, b, c, cgl),
  //        new NgsTestUtils.TestParams(READ_1,
  //            TEMP_1, EXP_SHORT1TOPN),
  //            new NgsTestUtils.NgsFilterPartlyParams(OutputFilter.TOPN, 5));
  //  }

  private static final String EXP_LONG2NONE = ""
          + "template\tR\t0\t-3\t0\t2" + LS
          + "template\tR\t0\t0\t0\t2" + LS
          + "template\tR\t0\t1\t0\t2" + LS
          + "template\tR\t0\t7\t0\t2" + LS
          + "template\tF\t0\t6\t0\t2" + LS
          + "template\tF\t0\t9\t0\t2" + LS
          + "template\tF\t0\t13\t0\t2" + LS
          + "template\tR\t0\t8\t0\t2" + LS
          + "template\tR\t0\t13\t0\t2" + LS;

  private static final String READ_2 = ""
      + ">read_1" + LS + "TGCTGT" + LS;
  //private static final String READ_2_REVERSE = ""
  //        + ">read_1" + LS + "ACAGCA" + LS;

  private static final String TEMP_2 =
      ">template" + LS
      + "CAGGCAACTGCCACCTT" + LS;

  public void testLong2None() throws Exception {
    final int w = 2, a = 0, b = 0, c = 1;
    final boolean cgl = false;
    final int stepSize = 2;
    final boolean useLongReads = true;
    checkLong(new NgsMaskParamsGeneral(w, a, b, c, cgl),
        new NgsTestUtils.TestParams(useLongReads,
            stepSize, READ_2, TEMP_2, EXP_LONG2NONE, 6),
            new NgsTestUtils.NgsFilterPartlyParams(OutputFilter.NONE, 0), false);
  }


  //  public void testLong2Topn() {
  //    final int r = 6, w = 2, a = 0, b = 0, c = 1;
  //    final boolean cgl = false;
  //    final int stepSize = 2;
  //    final boolean useLongReads = true;
  //    checkLong(new NgsMaskParamsGeneral(w, a, b, c, cgl),
  //        new NgsTestUtils.TestParams(useLongReads,
  //            stepSize, READ_2, TEMP_2, EXP_LONG2NONE),
  //            new NgsTestUtils.NgsFilterPartlyParams(OutputFilter.TOPN, 5), false);
  //  }

  private static final String SEQ_DNA_A = ">x" + LS + "actg" + LS;

  //similar to NgsTest.testA1()
  public void testA1Step4() throws Exception {
    final LogRecord lr = new LogRecord();
    Diagnostic.setLogStream(lr);
    final int stepSize = 4;
    final boolean useLongReads = true;
    checkLong(new NgsMaskParamsGeneral(4, 0, 0, 1, false),
        new NgsTestUtils.TestParams(useLongReads,
            stepSize, SEQ_DNA_A, SEQ_DNA_A, "x\tF\t0\t1\t0\t0" + LS, 4),
            new NgsTestUtils.NgsFilterPartlyParams(OutputFilter.NONE, 5), false);
    assertTrue(lr.toString().contains("GappedOutput score function calls: "));
    assertTrue(lr.toString().contains("GappedOutput scoreMax function calls: "));
    Diagnostic.setLogStream();
  }


  private static final String SEQ_DNA_ODD_S = "" + ">bar" + LS  + "acgt" + LS + ">baaa" + LS + "aaaa";

  private static final String SEQ_DNA_ODD_Q = "" + ">fooo" + LS + "acgt" + LS + ">fiii" + LS + "aaaa";

  private static final String EXPECTED_ODD = "fooo\tF\t0\t1\t0\t0" + LS + "fooo\tR\t0\t1\t0\t0" + LS + "fiii\tF\t1\t1\t0\t0" + LS;

  // only scores are different
  private static final String EXP_ODDLONG = "fooo\tF\t0\t1\t0\t0" + LS + "fooo\tR\t0\t1\t0\t0" + LS + "fiii\tF\t1\t1\t0\t0" + LS;

  // like NgsTest.testOddg
  public void testOddgLong() throws Exception {
    final int w = 4, a = 0, b = 0, c = 1;
    final boolean cgl = false;
    final int stepSize = 2;
    final boolean useLongReads = true;
    checkLong(new NgsMaskParamsGeneral(w, a, b, c, cgl),
        new NgsTestUtils.TestParams(useLongReads,
            stepSize, SEQ_DNA_ODD_S, SEQ_DNA_ODD_Q, EXP_ODDLONG, 8),
            new NgsTestUtils.NgsFilterPartlyParams(OutputFilter.NONE, 0), false);
  }

  // like NgsTest.testOddg
  public void testOddgShortNone() throws Exception {
    final int w = 4, a = 0, b = 0, c = 1;
    final boolean cgl = false;
    checkLong(new NgsMaskParamsGeneral(w, a, b, c, cgl),
        new NgsTestUtils.TestParams(SEQ_DNA_ODD_S,
            SEQ_DNA_ODD_Q, EXPECTED_ODD, 8),
            new NgsTestUtils.NgsFilterPartlyParams(OutputFilter.NONE, 0), false);
  }

  // like NgsTest.testOddg
  /**
   * @throws Exception
   */
  public void testOddgLongNone() throws Exception {
    final int w = 4, a = 0, b = 0, c = 1;
    final boolean cgl = false;
    final int stepSize = 2;
    final boolean useLongReads = true;
    checkLong(new NgsMaskParamsGeneral(w, a, b, c, cgl),
        new NgsTestUtils.TestParams(useLongReads,
            stepSize, SEQ_DNA_ODD_S, SEQ_DNA_ODD_Q, EXP_ODDLONG, 8),
            new NgsTestUtils.NgsFilterPartlyParams(OutputFilter.NONE, 0), false);
    //NgsTestUtils.check(new NgsMaskParamsGeneral(4, 4, 0, 0, 1, false), SEQ_DNA_ODD_S, SEQ_DNA_ODD_Q, "fooo\tF\t0\t1\t0\t0" + LS + "fooo\tR\t0\t1\t0\t0" + LS + "fiii\tF\t1\t1\t0\t0" + LS);
  }

  // like above with topn
  //  public void testOddgShortTopn() {
  //    final int r = 4, w = 4, a = 0, b = 0, c = 1;
  //    final boolean cgl = false;
  //    checkLong(new NgsMaskParamsGeneral(w, a, b, c, cgl),
  //        new NgsTestUtils.TestParams(SEQ_DNA_ODD_S,
  //            SEQ_DNA_ODD_Q, EXPECTED_ODD),
  //            new NgsTestUtils.NgsFilterPartlyParams(OutputFilter.TOPN, 5), false);
  //    //NgsTestUtils.check(new NgsMaskParamsGeneral(4, 4, 0, 0, 1, false), SEQ_DNA_ODD_S, SEQ_DNA_ODD_Q, "fooo\tF\t0\t1\t0\t0" + LS + "fooo\tR\t0\t1\t0\t0" + LS + "fiii\tF\t1\t1\t0\t0" + LS);
  //  }
  //
  //  // like above with topn
  //  public void testOddgLongTopn() {
  //    final int r = 4, w = 4, a = 0, b = 0, c = 1;
  //    final boolean cgl = false;
  //    final int stepSize = 2;
  //    final boolean useLongReads = true;
  //    checkLong(new NgsMaskParamsGeneral(w, a, b, c, cgl),
  //        new NgsTestUtils.TestParams(useLongReads,
  //            stepSize, SEQ_DNA_ODD_S, SEQ_DNA_ODD_Q, EXP_ODDLONG),
  //            new NgsTestUtils.NgsFilterPartlyParams(OutputFilter.TOPN, 5), false);
  //    //NgsTestUtils.check(new NgsMaskParamsGeneral(4, 4, 0, 0, 1, false), SEQ_DNA_ODD_S, SEQ_DNA_ODD_Q, "fooo\tF\t0\t1\t0\t0" + LS + "fooo\tR\t0\t1\t0\t0" + LS + "fiii\tF\t1\t1\t0\t0" + LS);
  //  }

  private static final String SEQ_DNA_ODD2_T = "" + ">bar" + LS  + "ggtacgtagg" + LS;

  private static final String SEQ_DNA_ODD3_T = "" + ">bar" + LS  + "tatacgtata" + LS;

  private static final String SEQ_DNA_ODD2_R = "" + ">fooo" + LS + "tacgta" + LS;

  private static final String EXP_ODD2 = ""
          + "bar\tF\t0\t-1\t0\t2" + LS
          + "bar\tR\t0\t-1\t0\t2" + LS
          + "bar\tF\t0\t7\t0\t2" + LS
          + "bar\tF\t0\t3\t0\t0" + LS
          + "bar\tR\t0\t3\t0\t0" + LS
          + "bar\tR\t0\t7\t0\t2" + LS;

    private static final String EXP_ODD2_W6 = ""
          + "bar\tF\t0\t3\t0\t0" + LS
          + "bar\tR\t0\t3\t0\t0" + LS;

  private static final String EXP_ODD3 = ""
          + "bar\tF\t0\t1\t0\t2" + LS
          + "bar\tF\t0\t-3\t0\t2" + LS
          + "bar\tR\t0\t-3\t0\t2" + LS
          + "bar\tR\t0\t1\t0\t2" + LS
          + "bar\tF\t0\t-1\t0\t2" + LS
          + "bar\tR\t0\t-1\t0\t2" + LS
          + "bar\tF\t0\t7\t0\t2" + LS
          + "bar\tF\t0\t3\t0\t0" + LS
          + "bar\tF\t0\t9\t0\t2" + LS
          + "bar\tF\t0\t5\t0\t2" + LS
          + "bar\tR\t0\t3\t0\t0" + LS
          + "bar\tR\t0\t7\t0\t2" + LS
          + "bar\tR\t0\t5\t0\t2" + LS
          + "bar\tR\t0\t9\t0\t2" + LS
      ;
  //
  //  public void testOdd2LongTopn() {
  //    final int r = 6, w = 6, a = 0, b = 0, c = 1;
  //    final boolean cgl = false;
  //    final int stepSize = 2;
  //    final boolean useLongReads = true;
  //    checkLong(new NgsMaskParamsGeneral(w, a, b, c, cgl),
  //        new NgsTestUtils.TestParams(useLongReads,
  //            stepSize, SEQ_DNA_ODD2_R, SEQ_DNA_ODD2_T, EXP_ODD2),
  //            new NgsTestUtils.NgsFilterPartlyParams(OutputFilter.TOPN, 5), false);
  //  }
  //
  //  public void testOdd2ShortTopn() {
  //    final int r = 6, w = 6, a = 0, b = 0, c = 1;
  //    final boolean cgl = false;
  //    NgsTestUtils.checkShort(new NgsMaskParamsGeneral(w, a, b, c, cgl),
  //        new NgsTestUtils.TestParams(SEQ_DNA_ODD2_R,
  //            SEQ_DNA_ODD2_T, EXP_ODD2),
  //            new NgsTestUtils.NgsFilterPartlyParams(OutputFilter.TOPN, 5));
  //  }

  public void testOdd2LongNone() throws Exception {
    final int w = 6, a = 0, b = 0, c = 1;
    final boolean cgl = false;
    final int stepSize = 2;
    final boolean useLongReads = true;
    checkLong(new NgsMaskParamsGeneral(w, a, b, c, cgl),
        new NgsTestUtils.TestParams(useLongReads,
            stepSize, SEQ_DNA_ODD2_R, SEQ_DNA_ODD2_T, EXP_ODD2_W6, 6),
            new NgsTestUtils.NgsFilterPartlyParams(OutputFilter.NONE, 0), false);
  }

  public void testOdd2LongNoneW2() throws Exception {
    final int w = 2, a = 0, b = 0, c = 1;
    final boolean cgl = false;
    final int stepSize = 2;
    final boolean useLongReads = true;
    checkLong(new NgsMaskParamsGeneral(w, a, b, c, cgl),
        new NgsTestUtils.TestParams(useLongReads,
            stepSize, SEQ_DNA_ODD2_R, SEQ_DNA_ODD2_T, EXP_ODD2, 6),
            new NgsTestUtils.NgsFilterPartlyParams(OutputFilter.NONE, 0), false);
  }

  public void testOdd2ShortNone() throws Exception {
    final int w = 6, a = 0, b = 0, c = 1;
    final boolean cgl = false;
    NgsTestUtils.checkShort(new NgsMaskParamsGeneral(w, a, b, c, cgl),
        new NgsTestUtils.TestParams(SEQ_DNA_ODD2_R,
            SEQ_DNA_ODD2_T, EXP_ODD2_W6, 6),
            new NgsTestUtils.NgsFilterPartlyParams(OutputFilter.NONE, 0));
  }

  public void testOdd3LongNone() throws Exception {
    final int w = 2, a = 0, b = 0, c = 1;
    final boolean cgl = false;
    final int stepSize = 2;
    final boolean useLongReads = true;
    checkLong(new NgsMaskParamsGeneral(w, a, b, c, cgl),
        new NgsTestUtils.TestParams(useLongReads,
            stepSize, SEQ_DNA_ODD2_R, SEQ_DNA_ODD3_T, EXP_ODD3, 6),
            new NgsTestUtils.NgsFilterPartlyParams(OutputFilter.NONE, 0), false);
  }

  static final String READS_BUG = ">read53:gi|48994873|gb|U00096.2|:4415055:S0:I1:D0" + LS + "ACCGTGCGTAATTTTTTATCACGGCTTaTACTTCAT" + LS;
  // 123123123123123123123123123123123123
  //"ACCGTGCGTAATTTTTTATCACGGCTTaTACTTCAT"
  //                            1  2  3  possible errors
  // 123456789123456789123456789123456789
  //                            1 possible error

  //                                                                                   aTAVTTCAT
  private static final String TEM_BUG = ">tbug" + LS + "cgttACCGTGCGTAATTTTTTATCACGGCTTTACTTCATcc" + LS;

  private static final String EXP_BUG_LONGS3 = "tbug\tF\t0\t5\t0\t2" + LS;

  private static final String EXP_BUG_LONGS9 = "tbug\tF\t0\t5\t0\t1" + LS;

  /** Like NgsTest.testBug36g */
  public void testBug36gLong() throws Exception {
    final int w = 18, a = 2, b = 1, c = 1;
    final boolean cgl = false;
    final int stepSize = 3;
    final boolean useLongReads = true;
    checkLong(new NgsMaskParamsGeneral(w, a, b, c, cgl),
        new NgsTestUtils.TestParams(useLongReads,
            stepSize, READS_BUG, TEM_BUG, EXP_BUG_LONGS3, 36),
            new NgsTestUtils.NgsFilterPartlyParams(OutputFilter.NONE, 5), false);
  }

  /** Like NgsTest.testBug36g */
  public void testBug36gLongStep9() throws Exception {
    final int w = 18, a = 2, b = 1, c = 1;
    final boolean cgl = false;
    final int stepSize = 9;
    final boolean useLongReads = true;
    checkLong(new NgsMaskParamsGeneral(w, a, b, c, cgl),
        new NgsTestUtils.TestParams(useLongReads,
            stepSize, READS_BUG, TEM_BUG, EXP_BUG_LONGS9, 36),
            new NgsTestUtils.NgsFilterPartlyParams(OutputFilter.NONE, 5), false);
  }

  // deletion
  private static final String READ_3 = ">read_1" + LS + "AAGG" + LS;

  private static final String TEMP_3 = ">template" + LS + "CCaTT" + LS;

  private static final String EXP_READ3 = "template\tR\t0\t1\t0\t2" + LS;
  //template  R 0 1 0 2

  public void testRead3Long() throws Exception {
    final int w = 2, a = 1, b = 1, c = 1;
    final boolean cgl = false;
    final int stepSize = 1;
    final boolean useLongReads = true;
    checkLong(new NgsMaskParamsGeneral(w, a, b, c, cgl),
        new NgsTestUtils.TestParams(useLongReads,
            stepSize, READ_3, TEMP_3, EXP_READ3, 4),
            new NgsTestUtils.NgsFilterPartlyParams(OutputFilter.NONE, 0), false);
    //    checkLong(new NgsMaskParamsGeneral(w, a, b, c, cgl),
    //        new NgsTestUtils.TestParams(useLongReads,
    //            stepSize, READ_3, TEMP_3, EXP_READ3),
    //            new NgsTestUtils.NgsFilterPartlyParams(OutputFilter.TOPN, 5), false);
  }

  // deletion
  private static final String READ_3F = ">read_1" + LS + "CCTT" + LS;

  private static final String TEMP_3F = ">template" + LS + "CCaTT" + LS;

  private static final String EXP_READ3F = "template\tF\t0\t1\t0\t2" + LS;
  //template  F 0 1 0 2

  // same result as reverse example
  public void testRead3FLong() throws Exception {
    final int w = 2, a = 1, b = 1, c = 1;
    final boolean cgl = false;
    final int stepSize = 1;
    final boolean useLongReads = true;
    checkLong(new NgsMaskParamsGeneral(w, a, b, c, cgl),
        new NgsTestUtils.TestParams(useLongReads,
            stepSize, READ_3F, TEMP_3F, EXP_READ3F, 4),
            new NgsTestUtils.NgsFilterPartlyParams(OutputFilter.NONE, 0), false);
    //    checkLong(new NgsMaskParamsGeneral(w, a, b, c, cgl),
    //        new NgsTestUtils.TestParams(useLongReads,
    //            stepSize, READ_3F, TEMP_3F, EXP_READ3F),
    //            new NgsTestUtils.NgsFilterPartlyParams(OutputFilter.TOPN, 5), false);
  }

  // insertion
  private static final String READ_4 = ">read_1" + LS + "AAccGG" + LS;
  //12345678901234567890123456789012345678901234567890123456789012

  private static final String TEMP_4 = ">template" + LS + "tttCCTTggg" + LS;

  private static final String EXP_READ4 = ""
          + "template\tR\t0\t-3\t0\t3" + LS
          + "template\tR\t0\t-2\t0\t3" + LS
          + "template\tF\t0\t2\t0\t6" + LS
          + "template\tF\t0\t5\t0\t3" + LS
          + "template\tR\t0\t2\t0\t3" + LS
          + "template\tR\t0\t4\t0\t6" + LS
          + "template\tR\t0\t7\t0\t4" + LS; //? My calculations say there should also be hits at F 4 and R 6, they don't seem to make it through gapped output

  //template  F 0 2 0 6
  //template  F 0 5 0 3
  //template  R 0 2 0 3
  //template  R 0 4 0 6

  // don't completely understand that one
  public void testRead4Long() throws Exception {
    final int w = 2, a = 2, b = 2, c = 1;
    final boolean cgl = false;
    final int stepSize = 1;
    final boolean useLongReads = true;
    final boolean ignoreOrder = false;
    checkLong(new NgsMaskParamsGeneral(w, a, b, c, cgl),
        new NgsTestUtils.TestParams(useLongReads,
            stepSize, READ_4, TEMP_4, EXP_READ4, 6),
            new NgsTestUtils.NgsFilterPartlyParams(OutputFilter.NONE, 0), ignoreOrder);
    //    ignoreOrder = true;
    //    checkLong(new NgsMaskParamsGeneral(w, a, b, c, cgl),
    //        new NgsTestUtils.TestParams(useLongReads,
    //            stepSize, READ_4, TEMP_4, EXP_READ4),
    //            new NgsTestUtils.NgsFilterPartlyParams(OutputFilter.TOPN, 5), ignoreOrder);
  }

  // insertion long
  private static final String READ_5 = ">read_1" + LS + "AAAAAAccGGGG" + LS;

  private static final String TEMP_5 = ">template" + LS + "AAAAAAGGGGtt" + LS;

  private static final String TEMP_5_1 = ">template" + LS + "ttAAAAAAGGGG" + LS;

  private static final String TEMP_5_2 = ">template" + LS + "CCCCTTTTTTtt" + LS;

  private static final String TEMP_5_3 = ">template" + LS + "ttCCCCTTTTTT" + LS;

  private static final String EXP_READ5 = ""
          + "template\tF\t0\t2\t0\t2" + LS
          + "template\tF\t0\t3\t0\t2" + LS
          + "template\tF\t0\t1\t0\t4" + LS
          + "template\tR\t0\t5\t0\t2" + LS;

  private static final String EXP_READ5_1 = ""
          + "template\tF\t0\t4\t0\t2" + LS
          + "template\tF\t0\t5\t0\t2" + LS
          + "template\tF\t0\t3\t0\t4" + LS;

  private static final String EXP_READ5_2 = ""
          + "template\tR\t0\t-3\t0\t2" + LS
          + "template\tR\t0\t-2\t0\t2" + LS
          + "template\tR\t0\t-1\t0\t2" + LS
          + "template\tR\t0\t0\t0\t2" + LS
          + "template\tR\t0\t1\t0\t1" + LS;

  private static final String EXP_READ5_3 = ""
          + "template\tR\t0\t-1\t0\t3" + LS
          + "template\tR\t0\t0\t0\t3" + LS
          + "template\tR\t0\t2\t0\t3" + LS
          + "template\tR\t0\t1\t0\t2" + LS
          + "template\tR\t0\t3\t0\t2" + LS;

  public void testRead5Long() throws Exception {
    final int w = 4, a = 2, b = 2, c = 1;
    final boolean cgl = false;
    final int stepSize = 4;
    final boolean useLongReads = true;
    checkLong(new NgsMaskParamsGeneral(w, a, b, c, cgl),
        new NgsTestUtils.TestParams(useLongReads,
            stepSize, READ_5, TEMP_5, EXP_READ5, 12),
            new NgsTestUtils.NgsFilterPartlyParams(OutputFilter.NONE, 0), false);
    //    checkLong(new NgsMaskParamsGeneral(w, a, b, c, cgl),
    //        new NgsTestUtils.TestParams(useLongReads,
    //            stepSize, READ_5, TEMP_5, EXP_READ5),
    //            new NgsTestUtils.NgsFilterPartlyParams(OutputFilter.TOPN, 5), false);
  }

  // stepsize 4 found nothing
  public void testRead5Long1() throws Exception {
    final int w = 4, a = 2, b = 2, c = 1;
    final boolean cgl = false;
    final int stepSize = 4;
    final boolean useLongReads = true;
    checkLong(new NgsMaskParamsGeneral(w, a, b, c, cgl),
        new NgsTestUtils.TestParams(useLongReads,
            stepSize, READ_5, TEMP_5_1, EXP_READ5_1, 12),
            new NgsTestUtils.NgsFilterPartlyParams(OutputFilter.NONE, 0), false);
    //    checkLong(new NgsMaskParamsGeneral(w, a, b, c, cgl),
    //        new NgsTestUtils.TestParams(useLongReads,
    //            stepSize, READ_5, TEMP_5_1, EXP_READ5_1),
    //            new NgsTestUtils.NgsFilterPartlyParams(OutputFilter.TOPN, 5), false);
  }
  public void testRead5Long2() throws Exception {
    final int w = 4, a = 2, b = 2, c = 1;
    final boolean cgl = false;
    final int stepSize = 4;
    final boolean useLongReads = true;
    checkLong(new NgsMaskParamsGeneral(w, a, b, c, cgl),
        new NgsTestUtils.TestParams(useLongReads,
            stepSize, READ_5, TEMP_5_2, EXP_READ5_2, 12),
            new NgsTestUtils.NgsFilterPartlyParams(OutputFilter.NONE, 0), false);
    //    checkLong(new NgsMaskParamsGeneral(w, a, b, c, cgl),
    //        new NgsTestUtils.TestParams(useLongReads,
    //            stepSize, READ_5, TEMP_5_2, EXP_READ5_2),
    //            new NgsTestUtils.NgsFilterPartlyParams(OutputFilter.TOPN, 5), false);
  }

  // no results with stepsize 4
  public void testRead5Long3() throws Exception {
    final int w = 4, a = 2, b = 2, c = 1;
    final boolean cgl = false;
    final int stepSize = 2;
    final boolean useLongReads = true;
    checkLong(new NgsMaskParamsGeneral(w, a, b, c, cgl),
        new NgsTestUtils.TestParams(useLongReads,
            stepSize, READ_5, TEMP_5_3, EXP_READ5_3, 12),
            new NgsTestUtils.NgsFilterPartlyParams(OutputFilter.NONE, 0), false);
    //    checkLong(new NgsMaskParamsGeneral(w, a, b, c, cgl),
    //        new NgsTestUtils.TestParams(useLongReads,
    //            stepSize, READ_5, TEMP_5_3, EXP_READ5_3),
    //            new NgsTestUtils.NgsFilterPartlyParams(OutputFilter.TOPN, 5), false);
  }

  // insertion really long
  //0AAAA
  //2AAAA
  //4AACC
  //6cccc
  //8cccc
  //10cccc
  //12cccG
  //14cGGG
  private static final String READ_6 = ">read_1" + LS +   "AAAAAAcccccccccGGGG" + LS;
  //                                                       01234567890123456789
  private static final String TEMP_6 = ">template" + LS + "AAAAAAGGGGccccccccc" + LS;
  //                                                         CCCCGGGGGGGGGTTTTTT
  //                                                       AAAAAACCCCCCCCCGGGG
  //                                                       gggggggggcccctttttt

  private static final String EXP_READ6 = ""
          + "template\tF\t0\t-1\t0\t5" + LS
          + "template\tF\t0\t0\t0\t5" + LS
          + "template\tF\t0\t1\t0\t4" + LS
          + "template\tF\t0\t2\t0\t4" + LS
          + "template\tF\t0\t3\t0\t4" + LS
          + "template\tF\t0\t4\t0\t5" + LS
          + "template\tF\t0\t9\t0\t5" + LS
          + "template\tF\t0\t7\t0\t5" + LS
          + "template\tF\t0\t5\t0\t4" + LS
          + "template\tF\t0\t10\t0\t5" + LS
          + "template\tF\t0\t8\t0\t5" + LS
          + "template\tF\t0\t6\t0\t4" + LS
          + "template\tR\t0\t-2\t0\t5" + LS
          + "template\tR\t0\t0\t0\t6" + LS
          + "template\tR\t0\t2\t0\t5" + LS;

 //private static final String EXP_READ5 = "template\tF\t0\t1\t0\t4" + LS;

  // same stepsize 2 as 4
  public void testRead6Long() throws Exception {
    final int w = 4, a = 9, b = 9, c = 1;
    final boolean cgl = false;
    final int stepSize = 2;
    final boolean useLongReads = true;
    checkLong(new NgsMaskParamsGeneral(w, a, b, c, cgl),
        new NgsTestUtils.TestParams(useLongReads,
            stepSize, READ_6, TEMP_6, EXP_READ6, 19),
            new NgsTestUtils.NgsFilterPartlyParams(OutputFilter.NONE, 0), false);
    //    checkLong(new NgsMaskParamsGeneral(w, a, b, c, cgl),
    //        new NgsTestUtils.TestParams(useLongReads,
    //            stepSize, READ_6, TEMP_6, EXP_READ5),
    //            new NgsTestUtils.NgsFilterPartlyParams(OutputFilter.TOPN, 5), false);
  }

  // insertion really long
  //0AAAA
  //2AAAA
  //4AACC
  //6cccc
  //8cccc
  //10cccc
  //12cccG
  //14cGGG                                                       0123456789012345678
  private static final String READ_REVERSE = ">read_1" + LS +   "gggggggggaaaatttttt" + LS;
  //                                                             aaaaaattttccccccccc
  private static final String TEMP_REVERSE = ">template" + LS + "AAAAAAGGGGccccccccc" + LS;
  //                                                             gggggggggcccctttttt

  private static final String EXP_REVERSE = ""
          + "template\tF\t0\t7\t0\t5" + LS
          + "template\tF\t0\t5\t0\t5" + LS
          + "template\tF\t0\t3\t0\t6" + LS
          + "template\tR\t0\t-4\t0\t5" + LS
          + "template\tR\t0\t-3\t0\t5" + LS
          + "template\tR\t0\t-2\t0\t4" + LS
          + "template\tR\t0\t-1\t0\t4" + LS
          + "template\tR\t0\t0\t0\t3" + LS
          + "template\tR\t0\t2\t0\t4" + LS
          + "template\tR\t0\t4\t0\t6" + LS
          + "template\tR\t0\t1\t0\t3" + LS
          + "template\tR\t0\t3\t0\t5" + LS
          + "template\tR\t0\t5\t0\t6" + LS;

  // same stepsize 2 as 4
  public void testRead6LongReverse() throws Exception {
    final int w = 4, a = 9, b = 9, c = 1;
    final boolean cgl = false;
    final int stepSize = 2;
    final boolean useLongReads = true;
    checkLong(new NgsMaskParamsGeneral(w, a, b, c, cgl),
        new NgsTestUtils.TestParams(useLongReads,
            stepSize, READ_REVERSE, TEMP_REVERSE, EXP_REVERSE, 19),
            new NgsTestUtils.NgsFilterPartlyParams(OutputFilter.NONE, 0), false);
    //    checkLong(new NgsMaskParamsGeneral(w, a, b, c, cgl),
    //        new NgsTestUtils.TestParams(useLongReads,
    //            stepSize, READ_REVERSE, TEMP_REVERSE, EXP_REVERSE),
    //            new NgsTestUtils.NgsFilterPartlyParams(OutputFilter.TOPN, 5), false);
  }

  // building reads together wrongly
  private static final String READ_7 = ""
      + ">read_1" + LS +   "AAAAGGGG" + LS
      + ">read_2" + LS +   "CCCCTTTT" + LS;

  private static final String TEMP_7 = ">template" + LS + "AAAAGGGGCCCCTTTT" + LS;

  private static final String EXP_READ7 = ""
      + "template\tF\t0\t1\t0\t0" + LS
      + "template\tF\t1\t9\t0\t0" + LS
      + "template\tR\t1\t1\t0\t0" + LS
      + "template\tR\t0\t9\t0\t0" + LS
      ;

  //template  F 0 1 0 0
  //template  F 1 9 0 0
  //template  R 0 9 0 0
  //template  R 1 1 0 0

  // same stepsize 2 as 4
  public void testRead7Long() throws Exception {
    final int w = 4, a = 4, b = 4, c = 1;
    final boolean cgl = false;
    final int stepSize = 2;
    final boolean useLongReads = true;
    checkLong(new NgsMaskParamsGeneral(w, a, b, c, cgl),
        new NgsTestUtils.TestParams(useLongReads,
            stepSize, READ_7, TEMP_7, EXP_READ7, 16),
            new NgsTestUtils.NgsFilterPartlyParams(OutputFilter.NONE, 0), false);
  }

  // test if building reads together wrongly?
  private static final String READ_NNN = ""
      + ">read_1" + LS +   "AAAAGGGG" + LS
      + ">read_2" + LS +   "CCCCNNNN" + LS;

  private static final String TEMP_NNN = ">template" + LS + "NNNNGGGGCCCCNNNN" + LS;

  private static final String EXP_READNNN = ""
      + "template\tF\t0\t1\t0\t2" + LS
      + "template\tF\t1\t9\t0\t2" + LS
      + "template\tR\t1\t1\t0\t2" + LS
      + "template\tR\t0\t9\t0\t2" + LS;


  // same stepsize 2 as 4
  public void testReadNNNLong() throws Exception {
    final int w = 4, a = 4, b = 4, c = 1;
    final boolean cgl = false;
    final int stepSize = 2;
    final boolean useLongReads = true;
    checkLong(new NgsMaskParamsGeneral(w, a, b, c, cgl),
        new NgsTestUtils.TestParams(useLongReads,
            stepSize, READ_NNN, TEMP_NNN, EXP_READNNN, 16),
            new NgsTestUtils.NgsFilterPartlyParams(OutputFilter.NONE, 0), false);
  }

  /// does not match with indel score 0 any more
  private static final String READ_NNN1 = ""
      + ">read_1" + LS +   "NGGGCCCC" + LS;

  private static final String TEMP_NNN1 = ">template" + LS + "NNNNGGGCCCCNNNN" + LS;

  private static final String EXP_READNNN1 = ""
      + "template\tF\t0\t4\t0\t1" + LS
      + "template\tR\t0\t4\t0\t2" + LS;

  // same stepsize 2 as 4
  public void testReadNNN1Long() throws Exception {
    final int w = 4, a = 4, b = 4, c = 1;
    final boolean cgl = false;
    final int stepSize = 1;
    final boolean useLongReads = true;
    boolean ignoreOrder = false;
    checkLong(new NgsMaskParamsGeneral(w, a, b, c, cgl),
        new NgsTestUtils.TestParams(useLongReads,
            stepSize, READ_NNN1, TEMP_NNN1, EXP_READNNN1, 8),
            new NgsTestUtils.NgsFilterPartlyParams(OutputFilter.NONE, 0), ignoreOrder);
    ignoreOrder = true;
    //    checkLong(new NgsMaskParamsGeneral(w, a, b, c, cgl),
    //        new NgsTestUtils.TestParams(useLongReads,
    //            stepSize, READ_NNN1, TEMP_NNN1, EXP_READNNN1),
    //            new NgsTestUtils.NgsFilterPartlyParams(OutputFilter.TOPN, 5), ignoreOrder);
  }

  //  private static final String READ_AMB = ""
  //    + ">read_1" + LS + "GCG" + LS;

  //  private static final String TEMP_AMB = ">template" + LS + "NNGCGNNNGCGNNNGCGNNNNNNGCGNNN" + LS;

  //  private static final String EXP_READAMB = ""
  //    + "#read 0 had 4 results with score-indel 0" + LS;

  // same stepsize 2 as 4
  //  public void testReadAmbLongTopn() {
  //    final int r = 3, w = 3, a = 0, b = 0, c = 1;
  //    final boolean cgl = false;
  //    final int stepSize = 1;
  //    final boolean useLongReads = true;
  //    final boolean ignoreOrder = false;
  //    checkLong(new NgsMaskParamsGeneral(w, a, b, c, cgl),
  //        new NgsTestUtils.TestParams(useLongReads,
  //            stepSize, READ_AMB, TEMP_AMB, EXP_READAMB),
  //            new NgsTestUtils.NgsFilterPartlyParams(OutputFilter.TOPN, 3), ignoreOrder);
  //  }


  void checkLong(final NgsMaskParams mask, final TestParams test, final NgsTestUtils.NgsFilterPartlyParams filterParams, final boolean containsOnly) throws Exception {
    NgsTestUtils.checkLong(mask, test, filterParams, containsOnly);
  }


}
