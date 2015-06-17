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
package com.rtg.alignment;

import static com.rtg.util.StringUtils.LS;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Locale;

import com.rtg.mode.DnaUtils;
import com.rtg.ngs.tempstage.BinaryTempFileRecord;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.PortableRandom;
import com.rtg.util.StringUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.variant.MachineErrorParams;
import com.rtg.variant.MachineErrorParamsBuilder;
import com.rtg.variant.realign.RealignParams;
import com.rtg.variant.realign.RealignParamsImplementation;
import com.rtg.variant.realign.ScoreMatrixCGTest;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 */
public class CgGotohEditDistanceTest extends TestCase {

  private RealignParams mParams;
  private CgGotohEditDistance mCgGotoh;
  private String mLog;
  private AlignmentResult mAlignment;
  private BinaryTempFileRecord mSamRecord;

  public static Test suite() {
    final TestSuite suite = new TestSuite();
    suite.addTestSuite(CgGotohEditDistanceTest.class);
    return suite;
  }

  @Override
  public void setUp() {
  }

  @Override
  public void tearDown() {
    mParams = null;
    mCgGotoh = null;
    mLog = null;
    mAlignment = null;
    mSamRecord = null;
  }

  protected EditDistance getEditDistanceInstance(int gapOpenPenalty, int gapExtendPenalty, int unknownsPenalty, String errorsFile) {
    try {
      mParams = new RealignParamsImplementation(new MachineErrorParamsBuilder().errors(errorsFile).create());
    } catch (final InvalidParamsException | IOException e) {
      throw new RuntimeException(e);
    }
    mCgGotoh = new CgGotohEditDistance(7, mParams, unknownsPenalty);
    return new RcEditDistance(mCgGotoh);
  }

  protected EditDistance getEditDistanceInstance(int gapOpenPenalty, int gapExtendPenalty, int unknownsPenalty) {
    return getEditDistanceInstance(gapOpenPenalty, gapExtendPenalty, unknownsPenalty, "cg_test_errors-080412");
  }

  protected void checkCG(String read, String template, String x1, String x2,
      String m, String xCigar, int score, boolean first, boolean rc, int start, boolean left) {
    checkCG(read, template, x1, x2, m, xCigar, score, first, rc, start, left, 0);
  }

  protected void checkCG(String read, String template, String x1, String x2,
      String m, String xCigar, int score, boolean first, boolean rc, int start, boolean left, int unknownsPenalty) {
    final byte[] s1 = DnaUtils.encodeString(read.replaceAll(" ", "").toLowerCase(Locale.ROOT));
    final byte[] s2 = DnaUtils.encodeString(template.replaceAll(" ", "").toLowerCase(Locale.ROOT));
    final ByteArrayOutputStream bos = new ByteArrayOutputStream();
    Diagnostic.setLogStream(new PrintStream(bos));
    final EditDistance f = getEditDistanceInstance(1, 1, unknownsPenalty);
    final int[] v = f.calculateEditDistance(s1, s1.length, s2, start, rc, Integer.MAX_VALUE, 7, left);
    f.logStats();
        //System.out.println(mCgGotoh.toString());
        //System.out.println("start: " + ActionsHelper.zeroBasedTemplateStart(v) + " diff:" + (start - ActionsHelper.zeroBasedTemplateStart(v))
         //   + " score: " + ActionsHelper.alignmentScore(v));
        //System.out.println("actions: " + ActionsHelper.toString(v));
    mAlignment = new AlignmentResult(s1, v, s2);
    mAlignment.setIdentifyingInfo(first, rc);

    final String cigar = mAlignment.getCigarString(rc, false);
    //     System.out.println("Ex:" + xCigar + "|" + cigar);
    final String ss = mAlignment.toString();
    //    System.err.println("alignment result: " + ss);
    mAlignment.setRemainingOutput(1, 0);
    mSamRecord = mAlignment.toRecord(false, null, 0, false, false);
    assertEquals(xCigar, new String(mSamRecord.getCigarString()));
    assertEquals(ss, score, ActionsHelper.alignmentScore(v));
    assertEquals(xCigar, cigar);

    assertEquals(x1.replaceAll("\\.", ""), mAlignment.readString());
    assertEquals(m, ActionsHelper.toString(v));

//    assertEquals(ss, x1 + "\t" + x2 + "\t" + m, mAlignment.tabularString());
    Diagnostic.setLogStream();
    mLog = bos.toString();
    //System.err.println("log: " + mLog);

    // now check that we get roughly the same alignment if we flip read+template onto the opposite arm.
    // (this is not always true, but the scores should always be pretty close)
    final byte[] s1R = DnaUtils.encodeString(reverse(read).replaceAll(" ", "").toLowerCase(Locale.ROOT));
    final byte[] s2R = DnaUtils.encodeString(reverse(template).replaceAll(" ", "").toLowerCase(Locale.ROOT));
    assertEquals(s1.length, s1R.length);
    assertEquals(s2.length, s2R.length);
    final EditDistance fR = getEditDistanceInstance(1, 1, unknownsPenalty); // because f.logStats() kills it!
    final int startR = s2.length - (start + s1.length) - 5;
    final int[] vR = fR.calculateEditDistance(s1R, s1R.length, s2R, startR, rc, Integer.MAX_VALUE, 7, !left);
        //System.out.println(mCgGotoh.toString());
        //System.out.println("actions: " + ActionsHelper.toString(vR));
        //System.err.println("startR: " + startR);
    assertEquals(ss, score, ActionsHelper.alignmentScore(vR));
    // the actions strings should be nearly the same, but can have a short reversed segment due to an ambiguous insert/delete position.
    final String actsLeft = ActionsHelper.toString(v);
    final String actsLeftRev = reverse(actsLeft);
    final String acts = ActionsHelper.toString(vR);
    assertEquals(actsLeftRev.length(), acts.length());
    assertEquals(actsLeftRev, acts);
  }

  protected void checkCG(String a, String b, String x1, String x2, String m, String xCigar, int score, boolean rc, boolean left) {
    checkCG(a, b, x1, x2, m, xCigar, score, true, rc, 0, left);
  }

  private String reverse(String str) {
    final StringBuilder sb = new StringBuilder();
    for (int i = str.length() - 1; i >= 0; i--) {
      sb.append(str.charAt(i));
    }
    return sb.toString();
  }

  // this is the example in com/rtg/variant/realign/scorematrixtestCG.xls
  static final String CG_SPREADSHEET = ""
    + "ScoreMatrix                      |G                0         |G                1         |G                2         |G                3         |G                4         |G                5         |G                6         |A                7         |T                8         |A                9         |A               10         |A               11         |A               12         |A               13         |A               14         |G               15         |G               16         |C               17         |G               18         |A               19         |C               20         |A               21         |T               22         |G               23         |C               24         |C               25         |A               26         |A               27         |T               28         |G               29         |T               30         |G               31         |T               32         |C               33         |G               34         |C               35         |C               36         |T               37         |T               38         |T               39         |T               40         |T               41         |C               42         |A               43         |A               44         |C               45         |T               46         |T               47         |T               48         |C               49         |C               50         |G               51         |A               52         |T               53         |" + LS
    + "[  0]           6.66e-02 5.53e-05|          6.66e-02 5.53e-05|          6.66e-02 5.53e-05|          6.66e-02 5.53e-05|          6.66e-02 5.53e-05|          6.66e-02 5.53e-05|          6.66e-02 5.53e-05|          6.66e-02 5.53e-05|          6.66e-02 5.53e-05|          6.66e-02 5.53e-05|          6.66e-02 5.53e-05|          6.66e-02 5.53e-05|          6.66e-02 5.53e-05|          6.66e-02 5.53e-05|          6.66e-02 5.53e-05|" + LS
    + "[  1]A                           |          4.23e-04 1.63e-05| 2.45e-07 4.23e-04 1.63e-05A 2.94e-07 4.23e-04 1.63e-05| 3.04e-07 4.23e-04 1.63e-05| 3.06e-07 4.23e-04 1.63e-05| 3.07e-07 4.23e-04 1.63e-05| 3.07e-07 4.23e-04 1.63e-05A 3.07e-07 6.52e-02 1.63e-05| 3.79e-05 4.23e-04 1.63e-05| 7.82e-06 6.52e-02 1.63e-05| 3.94e-05 6.52e-02 1.63e-05| 4.57e-05 6.52e-02 1.63e-05A 4.70e-05 6.52e-02 1.63e-05| 4.72e-05 6.52e-02 1.63e-05| 4.73e-05 6.52e-02         |" + LS
    + "[  2]T                           |                           |          2.77e-06 8.22e-07T 1.61e-09 2.77e-06 8.22e-07| 1.93e-09 2.77e-06 8.22e-07| 1.99e-09 2.77e-06 8.22e-07| 2.01e-09 2.77e-06 8.22e-07| 2.01e-09 2.77e-06 8.22e-07T 2.01e-09 2.77e-06 1.43e-05| 2.01e-09 6.38e-02 8.22e-07| 3.70e-05 2.96e-06 1.43e-05| 7.40e-06 4.14e-04 1.43e-05| 1.72e-06 4.14e-04 1.43e-05T 5.84e-07 4.14e-04 1.43e-05| 3.57e-07 4.14e-04 1.43e-05| 3.12e-07 4.14e-04 1.35e-05| 3.02e-07 4.14e-04         |" + LS
    + "[  3]A                           |                           |                           A          2.19e-08 3.76e-08| 1.27e-11 2.19e-08 3.76e-08| 1.52e-11 2.19e-08 3.76e-08| 1.57e-11 2.19e-08 3.76e-08| 1.58e-11 2.19e-08 3.76e-08A 1.59e-11 3.37e-06 6.42e-07| 1.96e-09 9.19e-08 1.33e-05| 4.45e-10 6.24e-02 6.42e-07| 3.62e-05 4.33e-05 7.28e-07| 7.26e-06 4.22e-04 7.28e-07A 1.70e-06 4.18e-04 7.28e-07| 5.82e-07 4.17e-04 7.28e-07| 3.58e-07 4.17e-04 6.95e-07| 3.13e-07 2.70e-06 8.59e-08| 6.42e-08 2.63e-06         A" + LS
    + "[  4]A                           |                           |                           A                           |          3.34e-10 1.69e-09| 1.94e-13 3.35e-10 1.69e-09| 2.33e-13 3.35e-10 1.69e-09| 2.41e-13 3.35e-10 1.69e-09A 2.42e-13 5.16e-08 2.96e-08| 3.00e-11 2.47e-08 5.97e-07| 2.03e-11 1.08e-05 1.30e-05| 6.24e-09 6.10e-02 4.17e-08| 3.54e-05 7.13e-05 1.20e-07A 7.12e-06 4.19e-04 1.19e-07| 1.67e-06 4.11e-04 1.19e-07| 5.72e-07 4.09e-04 1.18e-07| 3.51e-07 2.65e-06 4.43e-09| 7.18e-08 1.92e-08 5.45e-10A 1.44e-08 1.70e-08         |" + LS
    + "[  5]A                           |                           |                           A                           |                           |          1.10e-11 7.63e-11| 6.35e-15 1.10e-11 7.63e-11| 7.62e-15 1.10e-11 7.63e-11A 7.88e-15 1.69e-09 1.34e-09| 9.81e-13 4.82e-10 2.69e-08| 4.75e-13 5.04e-07 5.86e-07| 2.92e-10 2.09e-05 1.27e-05| 1.22e-08 5.97e-02 2.02e-08A 3.46e-05 9.76e-05 9.23e-08| 6.98e-06 4.16e-04 9.06e-08| 1.64e-06 4.03e-04 9.01e-08| 5.61e-07 2.60e-06 7.49e-10| 1.14e-07 1.86e-08 2.85e-11A 2.28e-08 4.90e-10 3.53e-12| 4.55e-09 1.81e-10         |" + LS
    + "      ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" + LS
    + "[  6]A                           |                           |                           A                           |          5.10e-13 1.81e-14| 2.96e-16 5.54e-13 6.70e-14| 3.81e-16 2.05e-12 9.84e-13| 1.26e-15 3.01e-11 2.29e-11A 1.77e-14 1.08e-07 7.48e-10| 6.27e-11 2.29e-08 9.97e-07| 2.58e-11 4.70e-03 1.04e-05| 2.73e-06 4.91e-02 1.01e-06| 2.90e-05 4.81e-03 8.08e-08A 8.59e-06 3.88e-04 7.72e-08| 1.94e-06 3.65e-04 7.14e-09| 6.00e-07 3.42e-05 4.64e-11| 1.40e-07 2.14e-09 3.99e-13| 2.80e-08 1.57e-10 3.97e-14A 5.60e-09 2.99e-11 3.00e-15|" + LS
    + "[  7]A                           |                           |                           A                           |                           |          3.33e-15 3.13e-15| 1.93e-18 3.87e-15 4.47e-14| 2.63e-18 1.81e-14 1.04e-12A 1.10e-17 4.79e-11 5.61e-11| 2.78e-14 6.90e-10 4.49e-08| 4.06e-13 8.23e-07 1.44e-06| 4.77e-10 4.60e-03 1.02e-05| 2.67e-06 4.80e-02 1.00e-06A 2.84e-05 4.73e-03 8.39e-08| 8.41e-06 3.86e-04 7.62e-08| 1.91e-06 3.59e-04 7.09e-09| 5.90e-07 2.20e-07 4.62e-13| 1.18e-07 7.25e-10 3.45e-14A 2.36e-08 1.43e-10 6.34e-15| 4.72e-09 2.86e-11         |" + LS
    + "[  8]G                           |                           |                           G                           |                           |                           |          5.77e-15 2.01e-15| 3.35e-18 3.97e-14 4.67e-14G 2.37e-17 5.52e-15 2.53e-12| 7.94e-18 5.96e-13 2.02e-09| 3.47e-16 2.38e-10 6.51e-08| 1.38e-13 1.27e-08 1.42e-06| 7.42e-12 2.93e-05 1.00e-05G 1.70e-08 3.05e-04 9.84e-07| 1.80e-07 3.01e-05 8.35e-08| 5.35e-08 2.49e-06 7.48e-08| 1.21e-08 3.53e-04 4.56e-11| 2.07e-07 6.77e-07 1.52e-13G 4.18e-08 6.05e-10 3.00e-14| 8.36e-09 1.86e-08 5.94e-15| 1.68e-09 2.42e-11         |" + LS
    + "[  9]G                           |                           |                           G                           |                           |                           |                           |          7.26e-15 2.11e-15G 4.21e-18 4.95e-16 1.14e-13| 1.13e-18 1.32e-14 9.09e-11| 7.91e-18 1.05e-11 2.93e-09| 6.11e-15 3.41e-10 6.37e-08| 1.99e-13 7.46e-09 4.56e-07G 4.37e-12 2.38e-07 1.07e-07| 1.39e-10 1.94e-06 1.00e-08| 1.15e-09 1.93e-07 3.88e-09| 3.42e-10 2.54e-06 7.32e-08| 1.54e-09 3.45e-04 1.40e-10G 2.00e-07 5.35e-09 1.27e-13| 4.01e-08 3.33e-08 3.87e-12| 8.04e-09 1.61e-10 5.02e-15| 1.61e-09 8.71e-12         |" + LS
    + "[ 10]C                           |                           |                           C                           |                           |                           |                           |                           C          5.71e-17 5.13e-15| 3.31e-20 5.97e-16 4.09e-12| 3.53e-19 4.73e-13 1.32e-10| 2.75e-16 1.53e-11 2.87e-09| 8.95e-15 3.34e-10 2.05e-08C 1.96e-13 2.42e-09 4.89e-09| 1.45e-12 2.07e-09 8.53e-10| 1.49e-12 1.23e-08 2.15e-10| 7.46e-12 1.25e-09 3.82e-09| 2.22e-12 1.65e-08 7.16e-08C 1.00e-11 3.37e-04 1.12e-12| 1.96e-07 1.05e-09 7.09e-12| 3.91e-08 4.15e-10 3.36e-14| 7.83e-09 6.45e-09 1.81e-15| 1.57e-09 8.23e-12         C" + LS
    + "[ 11]G                           |                           |                           G                           |                           |                           |                           |                           G                           |          2.71e-17 1.84e-13| 1.57e-20 2.13e-14 5.93e-12| 1.24e-17 6.90e-13 1.29e-10| 4.03e-16 1.50e-11 9.24e-10G 8.80e-15 1.09e-10 2.20e-10| 6.50e-14 4.08e-11 3.88e-11| 3.67e-14 1.76e-11 1.22e-11| 1.75e-14 1.23e-08 1.72e-10| 7.11e-12 4.29e-09 3.22e-09G 3.91e-12 4.78e-10 7.00e-08| 1.06e-12 3.30e-04 5.38e-13| 1.91e-07 1.00e-09 8.77e-14| 3.83e-08 2.02e-10 1.34e-12| 7.66e-09 8.08e-11 1.71e-15G 1.53e-09 8.03e-12         |" + LS
    + "[ 12]A                           |                           |                           A                           |                           |                           |                           |                           A                           |                           |          1.48e-13 2.67e-13| 8.57e-17 4.79e-12 5.81e-12| 2.79e-15 1.04e-10 4.16e-11A 6.11e-14 7.57e-10 9.94e-12| 4.51e-13 2.84e-10 1.75e-12| 2.55e-13 7.12e-11 5.54e-13| 9.22e-14 1.75e-13 1.03e-11| 1.85e-14 7.86e-11 1.46e-10A 4.93e-14 4.41e-11 3.15e-09| 3.54e-14 3.68e-10 6.85e-08| 2.20e-13 3.23e-04 2.12e-13| 1.87e-07 9.79e-10 1.02e-13| 3.74e-08 3.02e-08 1.68e-14A 7.51e-09 3.94e-11 1.67e-15| 1.50e-09 7.84e-12         |" + LS
    + "[ 13]C                           |                           |                           C                           |                           |                           |                           |                           C                           |                           |                           |          2.33e-15 2.62e-13| 1.35e-18 6.06e-14 1.89e-12C 3.54e-17 8.79e-13 6.04e-13| 5.17e-16 4.85e-12 1.38e-13| 2.92e-15 1.81e-12 3.97e-14| 1.63e-15 4.56e-13 4.63e-13| 5.91e-16 5.52e-14 6.59e-12C 1.50e-16 1.94e-10 1.42e-10| 1.13e-13 1.67e-11 3.08e-09| 3.22e-14 3.59e-10 6.70e-08| 2.15e-13 3.16e-04 2.08e-13| 1.83e-07 9.58e-10 6.27e-12C 3.66e-08 3.82e-10 8.26e-15| 7.33e-09 3.84e-11 1.63e-15| 1.47e-09 1.18e-09         |" + LS
    + "[ 14]A                           |                           |                           A                           |                           |                           |                           |                           A                           |                           |                           |                           |          2.13e-13 8.52e-14A 1.23e-16 1.58e-12 2.74e-14| 9.41e-16 1.34e-12 7.21e-15| 9.68e-16 4.86e-12 2.16e-15| 3.01e-15 1.17e-14 2.09e-14| 6.09e-16 5.31e-15 2.96e-13A 1.25e-16 3.47e-14 6.42e-12| 4.51e-17 1.97e-12 1.39e-10| 1.15e-15 2.49e-09 3.01e-09| 1.45e-12 3.51e-10 6.55e-08| 4.93e-13 3.09e-04 4.81e-13A 1.79e-07 9.37e-10 7.96e-14| 3.58e-08 1.89e-10 8.04e-15| 7.17e-09 3.75e-11 2.46e-13| 1.43e-09 1.50e-11         |" + LS
    + "[ 15]T                           |                           |                           T                           |                           |                           |                           |                           T                           |                           |                           |                           |                           T          1.80e-15 1.56e-15| 1.04e-18 1.02e-14 6.03e-16| 6.11e-18 8.58e-15 1.11e-15| 6.19e-18 3.09e-14 9.44e-16| 1.91e-17 1.99e-16 1.33e-14T 3.94e-18 1.58e-15 2.89e-13| 1.71e-18 3.37e-14 6.24e-12| 1.99e-17 7.35e-13 1.36e-10| 4.30e-16 3.15e-11 2.95e-09| 1.84e-14 3.44e-10 6.41e-08T 2.03e-13 3.02e-04 1.98e-13| 1.75e-07 9.17e-10 3.95e-14| 3.50e-08 1.83e-10 1.88e-14| 7.01e-09 3.67e-11 3.10e-15| 1.40e-09 7.38e-12         T" + LS
    + "      --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" + LS
    + "[ 16]G                           |                           |                           G                           |                           |                           |                           |                           G                           |                           |                           |                           |                           G                           |          5.75e-18 1.05e-18| 3.33e-21 3.20e-17 2.03e-18| 1.92e-20 9.56e-15 3.27e-18| 5.55e-18 1.54e-14 5.14e-18G 1.01e-17 1.57e-16 1.87e-17| 2.10e-18 8.82e-14 3.90e-16| 5.15e-17 1.19e-14 8.50e-15| 1.72e-17 2.60e-13 1.85e-13| 1.54e-16 5.66e-12 4.01e-12G 3.32e-15 1.23e-10 1.69e-08| 7.17e-14 7.98e-05 3.98e-08| 4.63e-08 1.22e-06 5.96e-09| 9.97e-09 1.83e-07 4.43e-14| 2.10e-09 2.09e-10 8.86e-15G 4.20e-10 4.17e-11 1.70e-15|" + LS
    + "[ 17]C                           |                           |                           C                           |                           |                           |                           |                           C                           |                           |                           |                           |                           C                           |                           |          4.19e-20 9.79e-20| 2.43e-23 2.13e-19 2.13e-18| 1.29e-22 6.07e-17 3.43e-18C 3.52e-20 1.51e-14 8.74e-19| 8.76e-18 1.15e-18 3.58e-17| 1.75e-18 5.61e-16 3.85e-16| 6.76e-19 1.85e-14 8.38e-15| 1.09e-17 2.62e-15 1.81e-13C 3.70e-18 5.68e-14 7.62e-10| 3.37e-17 8.90e-11 1.84e-08| 5.16e-14 7.81e-05 5.21e-10| 4.53e-08 1.23e-06 3.79e-11| 9.78e-09 1.21e-09 4.37e-14C 1.96e-09 1.20e-11 8.74e-15| 3.91e-10 2.40e-12         |" + LS
    + "[ 18]C                           |                           |                           C                           |                           |                           |                           |                           C                           |                           |                           |                           |                           C                           |                           |                           |          7.76e-22 9.59e-20| 4.50e-25 1.25e-20 1.67e-19C 7.32e-24 6.21e-17 3.17e-18| 3.60e-20 9.58e-17 1.61e-18| 6.27e-20 2.39e-19 1.74e-17| 1.27e-20 8.60e-16 3.81e-16| 5.01e-19 1.61e-16 8.16e-15C 1.94e-19 9.62e-16 3.43e-11| 5.97e-19 3.97e-12 8.26e-10| 2.30e-15 1.48e-08 1.62e-08| 8.60e-12 7.64e-05 2.57e-10| 4.43e-08 8.05e-09 2.53e-13C 8.87e-09 5.74e-11 2.88e-15| 1.77e-09 1.00e-11 4.98e-16| 3.55e-10 2.00e-12         |" + LS
    + "[ 19]A                           |                           |                           A                           |                           |                           |                           |                           A                           |                           |                           |                           |                           A                           |                           |                           |                           |          5.05e-22 7.52e-21A 2.93e-25 9.49e-22 1.56e-19| 6.09e-25 4.11e-19 9.25e-20| 2.38e-22 9.50e-17 7.85e-19| 5.51e-20 9.27e-20 1.73e-17| 1.11e-20 1.15e-15 3.67e-16A 6.68e-19 4.36e-17 1.54e-12| 1.59e-19 1.79e-13 3.72e-11| 1.04e-16 4.33e-12 7.34e-10| 2.53e-15 1.79e-10 1.59e-08| 1.04e-13 7.47e-05 1.68e-12A 4.34e-08 4.26e-08 1.20e-14| 8.69e-09 4.54e-11 2.10e-15| 1.74e-09 9.08e-12 4.16e-16| 3.48e-10 1.82e-12         |" + LS
    + "[ 20]A                           |                           |                           A                           |                           |                           |                           |                           A                           |                           |                           |                           |                           A                           |                           |                           |                           |                           A          4.24e-23 7.00e-21| 2.46e-26 8.17e-22 4.25e-21| 4.79e-25 4.76e-19 5.50e-20| 2.76e-22 6.07e-19 7.80e-19| 4.07e-22 1.41e-17 1.68e-17A 8.23e-21 9.19e-18 6.94e-14| 6.98e-21 8.04e-15 1.67e-12| 4.67e-18 1.95e-13 3.30e-11| 1.14e-16 3.85e-12 7.14e-10| 2.26e-15 1.29e-08 1.55e-08A 7.49e-12 7.31e-05 8.84e-12| 4.24e-08 4.91e-10 9.53e-15| 8.48e-09 4.45e-11 1.90e-15| 1.70e-09 8.90e-12 3.77e-16| 3.39e-10 1.78e-12         A" + LS
    + "[ 21]T                           |                           |                           T                           |                           |                           |                           |                           T                           |                           |                           |                           |                           T                           |                           |                           |                           |                           T                           |          3.67e-23 1.91e-22| 2.13e-26 2.73e-23 2.58e-21| 2.01e-26 3.31e-21 3.52e-20| 1.92e-24 7.92e-21 7.58e-19T 4.98e-24 2.72e-17 3.13e-15| 1.58e-20 3.62e-16 7.53e-14| 2.13e-19 8.77e-15 1.49e-12| 5.13e-18 1.73e-13 3.21e-11| 1.02e-16 3.75e-12 7.01e-10T 2.19e-15 1.63e-10 1.52e-08| 9.49e-14 7.15e-05 1.02e-13| 4.15e-08 2.19e-10 9.32e-15| 8.30e-09 6.69e-09 1.86e-15| 1.66e-09 8.68e-12 3.69e-16T 3.33e-10 2.68e-10         |" + LS
    + "[ 22]G                           |                           |                           G                           |                           |                           |                           |                           G                           |                           |                           |                           |                           G                           |                           |                           |                           |                           G                           |                           |          1.23e-24 1.16e-22| 7.13e-28 1.36e-23 1.59e-21| 8.03e-27 2.05e-22 3.41e-20G 1.20e-25 4.00e-21 1.41e-16| 2.34e-24 2.54e-15 3.39e-15| 1.47e-18 3.95e-16 6.69e-14| 5.23e-19 7.80e-15 1.45e-12| 4.63e-18 1.69e-13 3.15e-11G 9.87e-17 3.67e-12 6.83e-10| 2.15e-15 8.01e-11 1.48e-08| 4.69e-14 6.99e-05 4.58e-14| 4.06e-08 2.12e-10 1.39e-12| 8.11e-09 1.30e-08 1.82e-15G 1.63e-09 8.51e-12 5.55e-14| 3.26e-10 3.39e-12         |" + LS
    + "[ 23]T                           |                           |                           T                           |                           |                           |                           |                           T                           |                           |                           |                           |                           T                           |                           |                           |                           |                           T                           |                           |                           |          6.12e-25 7.14e-23| 3.55e-28 8.35e-24 1.53e-21T 4.91e-27 2.76e-20 6.33e-18| 1.60e-23 7.33e-19 1.53e-16| 4.28e-22 3.37e-17 3.01e-15| 1.97e-20 3.51e-16 6.51e-14| 2.07e-19 7.58e-15 1.42e-12T 4.44e-18 1.65e-13 3.07e-11| 9.68e-17 5.52e-10 6.68e-10| 3.20e-13 7.78e-11 1.45e-08| 1.09e-13 6.84e-05 1.07e-13| 3.97e-08 2.08e-10 2.71e-12T 7.94e-09 1.91e-08 4.26e-15| 1.60e-09 8.34e-12 7.03e-16| 3.20e-10 1.68e-12         |" + LS
    + "[ 24]G                           |                           |                           G                           |                           |                           |                           |                           G                           |                           |                           |                           |                           G                           |                           |                           |                           |                           G                           |                           |                           |                           |          3.76e-25 6.91e-23G 2.18e-28 8.05e-24 2.85e-19| 4.71e-27 5.11e-18 6.88e-18| 2.96e-21 8.02e-19 1.35e-16| 1.06e-21 1.59e-17 2.93e-15| 9.43e-21 3.41e-16 6.38e-14G 2.00e-19 7.44e-15 1.38e-12| 4.36e-18 1.61e-13 3.02e-11| 9.43e-17 1.08e-09 6.53e-10| 6.24e-13 7.61e-11 1.42e-08| 1.69e-13 6.69e-05 1.65e-13G 3.88e-08 2.03e-10 3.97e-12| 7.76e-09 1.62e-10 1.76e-15| 1.55e-09 1.26e-09 3.48e-16| 3.11e-10 1.64e-12         |" + LS
    + "[ 25]T                           |                           |                           T                           |                           |                           |                           |                           T                           |                           |                           |                           |                           T                           |                           |                           |                           |                           T                           |                           |                           |                           |                           T          5.58e-23 1.28e-20| 3.24e-26 1.48e-21 3.11e-19| 8.67e-25 6.83e-20 6.09e-18| 3.98e-23 7.11e-19 1.32e-16| 4.20e-22 1.54e-17 2.87e-15T 8.99e-21 3.35e-16 6.22e-14| 1.96e-19 1.12e-12 1.36e-12| 6.48e-16 1.58e-13 2.96e-11| 2.21e-16 1.58e-09 6.39e-10| 9.15e-13 7.45e-11 1.39e-08T 2.26e-13 6.55e-05 2.21e-13| 3.80e-08 1.99e-10 3.36e-14| 7.59e-09 4.05e-11 2.62e-13| 1.52e-09 1.59e-11 3.39e-16| 3.04e-10 1.59e-12         T" + LS
    + "      ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" + LS
    + "[ 26]T                           |                           |                           T                           |                           |                           |                           |                           T                           |                           |                           |                           |                           T                           |                           |                           |                           |                           T                           |                           |                           |                           |                           T                           |                           |                           |                           |                           T                           |                           |          1.79e-21 1.16e-21| 1.04e-24 5.45e-18 2.48e-20| 3.16e-21 7.60e-19 5.41e-19T 1.07e-21 2.55e-15 1.17e-17| 1.48e-18 3.59e-16 4.63e-16| 5.04e-19 1.42e-14 5.60e-15| 8.32e-18 1.71e-13 4.14e-13| 1.01e-16 1.27e-11 2.64e-12T 7.37e-15 1.25e-08 1.22e-08| 7.22e-12 5.76e-05 9.51e-10| 3.34e-08 4.51e-06 4.08e-10| 9.30e-09 1.93e-06 4.80e-15| 2.98e-09 2.40e-09 7.82e-16T 5.97e-10 3.11e-12 1.22e-16| 1.19e-10 3.43e-13 9.92e-18|" + LS
    + "[ 27]T                           |                           |                           T                           |                           |                           |                           |                           T                           |                           |                           |                           |                           T                           |                           |                           |                           |                           T                           |                           |                           |                           |                           T                           |                           |                           |                           |                           T                           |                           |                           |          2.68e-21 2.25e-21| 1.55e-24 3.47e-20 2.45e-20T 2.04e-23 1.18e-18 1.06e-18| 6.89e-22 1.63e-17 2.09e-17| 9.57e-21 4.70e-18 2.55e-16| 4.64e-21 1.19e-16 1.87e-14| 7.00e-20 3.25e-15 1.21e-13T 1.90e-18 1.45e-11 5.53e-10| 8.42e-15 2.20e-08 1.20e-08| 1.28e-11 5.64e-05 9.54e-10| 3.27e-08 4.44e-06 4.00e-10| 9.11e-09 1.89e-06 4.99e-13T 2.92e-09 3.04e-11 6.51e-16| 5.84e-10 3.06e-12 7.17e-17| 1.17e-10 6.09e-13         |" + LS
    + "[ 28]C                           |                           |                           C                           |                           |                           |                           |                           C                           |                           |                           |                           |                           C                           |                           |                           |                           |                           C                           |                           |                           |                           |                           C                           |                           |                           |                           |                           C                           |                           |                           |                           |          2.87e-23 1.11e-21C 1.67e-26 3.48e-22 4.78e-20| 2.05e-25 2.00e-18 9.45e-19| 1.16e-21 2.12e-19 1.15e-17| 3.56e-22 2.09e-16 8.41e-16| 1.21e-19 1.51e-14 5.46e-15C 8.79e-18 6.53e-16 2.49e-11| 2.14e-18 2.97e-12 5.45e-10| 1.72e-15 2.02e-10 1.17e-08| 1.18e-13 3.58e-07 9.39e-10| 2.07e-10 2.83e-08 3.93e-10C 5.79e-11 1.86e-06 6.34e-15| 1.09e-09 1.50e-11 6.37e-16| 2.18e-10 2.99e-12 1.26e-16| 4.36e-11 9.22e-11         |" + LS
    + "[ 29]A                           |                           |                           A                           |                           |                           |                           |                           A                           |                           |                           |                           |                           A                           |                           |                           |                           |                           A                           |                           |                           |                           |                           A                           |                           |                           |                           |                           A                           |                           |                           |                           |                           A          5.97e-24 2.15e-21| 3.46e-27 2.51e-22 4.29e-20| 1.47e-25 1.76e-20 5.16e-19| 1.03e-23 6.11e-20 3.79e-17| 3.75e-23 5.71e-18 2.49e-16A 3.32e-21 1.24e-16 1.12e-12| 7.28e-20 1.30e-13 2.45e-11| 7.52e-17 2.86e-12 5.28e-10| 1.67e-15 6.25e-11 1.16e-10| 3.66e-14 2.27e-09 2.36e-11A 1.33e-12 1.83e-10 3.86e-10| 3.71e-13 1.82e-06 3.15e-15| 1.05e-09 8.69e-10 6.26e-16| 2.11e-10 1.13e-12 1.91e-14| 4.23e-11 8.06e-13         A" + LS
    + "[ 30]A                           |                           |                           A                           |                           |                           |                           |                           A                           |                           |                           |                           |                           A                           |                           |                           |                           |                           A                           |                           |                           |                           |                           A                           |                           |                           |                           |                           A                           |                           |                           |                           |                           A                           |          1.13e-23 1.93e-21| 6.53e-27 2.25e-22 2.32e-20| 1.32e-25 2.80e-21 1.70e-18| 1.65e-24 1.98e-19 1.12e-17A 1.15e-22 1.33e-18 5.04e-14| 7.97e-22 5.83e-15 1.10e-12| 3.38e-18 1.29e-13 2.38e-11| 7.52e-17 2.77e-12 5.25e-12| 1.62e-15 1.00e-12 1.53e-12A 9.06e-16 1.46e-11 1.74e-11| 8.62e-15 4.90e-10 3.77e-10| 2.86e-13 1.78e-06 1.80e-13| 1.03e-09 1.09e-11 1.09e-15| 2.06e-10 1.08e-12 1.67e-16A 4.13e-11 2.20e-13         |" + LS
    + "[ 31]C                           |                           |                           C                           |                           |                           |                           |                           C                           |                           |                           |                           |                           C                           |                           |                           |                           |                           C                           |                           |                           |                           |                           C                           |                           |                           |                           |                           C                           |                           |                           |                           |                           C                           |                           |          1.01e-23 1.05e-21| 5.88e-27 1.89e-20 7.67e-20| 1.09e-23 1.37e-18 5.04e-19C 7.98e-22 5.97e-20 2.27e-15| 1.94e-22 2.62e-16 4.96e-14| 1.52e-19 5.78e-15 1.07e-12| 3.39e-18 1.25e-13 2.37e-13| 7.30e-17 4.50e-14 6.92e-14C 4.07e-17 2.21e-12 7.86e-13| 1.29e-15 1.83e-13 1.71e-11| 3.65e-16 5.07e-12 3.69e-10| 3.02e-15 1.74e-06 2.31e-15| 1.01e-09 5.31e-12 2.32e-16C 2.02e-10 1.06e-12 4.57e-17| 4.04e-11 2.11e-13         |" + LS
    + "[ 32]T                           |                           |                           T                           |                           |                           |                           |                           T                           |                           |                           |                           |                           T                           |                           |                           |                           |                           T                           |                           |                           |                           |                           T                           |                           |                           |                           |                           T                           |                           |                           |                           |                           T                           |                           |                           |          5.51e-24 3.46e-21| 3.20e-27 5.19e-22 2.30e-20T 3.02e-25 1.75e-18 1.02e-16| 1.01e-21 1.82e-15 2.23e-15| 1.06e-18 4.01e-14 4.81e-14| 2.35e-17 8.65e-13 1.07e-14| 5.06e-16 3.12e-13 3.12e-15T 2.82e-16 6.46e-16 3.58e-14| 5.69e-17 1.81e-14 7.69e-13| 2.19e-17 9.02e-14 1.66e-11| 5.67e-17 1.96e-12 3.61e-10| 1.15e-15 1.70e-06 1.11e-15T 9.87e-10 7.96e-10 2.21e-16| 1.98e-10 1.59e-10 4.38e-17| 3.97e-11 2.07e-13         |" + LS
    + "[ 33]T                           |                           |                           T                           |                           |                           |                           |                           T                           |                           |                           |                           |                           T                           |                           |                           |                           |                           T                           |                           |                           |                           |                           T                           |                           |                           |                           |                           T                           |                           |                           |                           |                           T                           |                           |                           |                           |          1.80e-23 1.03e-21T 1.05e-26 1.90e-20 4.59e-18| 1.10e-23 8.36e-17 1.01e-16| 4.85e-20 3.57e-15 2.18e-15| 2.08e-18 7.79e-14 6.61e-16| 4.56e-17 8.55e-13 2.05e-16T 5.05e-16 2.00e-15 1.61e-15| 1.02e-16 1.92e-16 3.46e-14| 2.05e-17 4.12e-15 7.48e-13| 6.50e-18 8.71e-14 1.62e-11| 5.18e-17 2.92e-10 3.53e-10T 1.69e-13 1.66e-06 1.65e-13| 9.66e-10 1.55e-09 3.30e-14| 1.94e-10 2.02e-12 4.29e-17| 3.88e-11 2.03e-13         |" + LS
    + "[ 34]T                           |                           |                           T                           |                           |                           |                           |                           T                           |                           |                           |                           |                           T                           |                           |                           |                           |                           T                           |                           |                           |                           |                           T                           |                           |                           |                           |                           T                           |                           |                           |                           |                           T                           |                           |                           |                           |                           T          8.48e-22 2.07e-19| 4.92e-25 3.71e-18 4.56e-18| 2.15e-21 1.63e-16 9.86e-17| 9.49e-20 5.24e-15 4.59e-17| 3.06e-18 7.67e-14 1.87e-16T 4.51e-17 5.43e-15 7.30e-17| 1.22e-17 2.37e-17 1.56e-15| 2.45e-18 1.82e-16 3.36e-14| 5.95e-19 3.92e-15 7.31e-13| 2.39e-18 1.31e-11 1.60e-11T 7.62e-15 5.69e-10 3.45e-10| 3.32e-13 1.63e-06 3.24e-13| 9.45e-10 1.48e-11 4.20e-16| 1.89e-10 9.99e-13 4.21e-17| 3.78e-11 1.99e-13         T" + LS
    + "[ 35]C                           |                           |                           C                           |                           |                           |                           |                           C                           |                           |                           |                           |                           C                           |                           |                           |                           |                           C                           |                           |                           |                           |                           C                           |                           |                           |                           |                           C                           |                           |                           |                           |                           C                           |                           |                           |                           |                           C                           |          1.08e-21 2.06e-19|          4.73e-20 4.47e-18|          1.55e-18 3.15e-18|          3.35e-17 2.43e-17C          7.52e-14 4.41e-18|          3.50e-17 7.01e-17|          8.33e-18 1.51e-15|          2.72e-14 3.29e-14|          3.83e-15 7.21e-13C          1.66e-13 1.57e-11|          5.41e-12 3.38e-10|          1.59e-06 3.08e-15|          7.55e-10 2.09e-16|          9.67e-13 4.12e-17C          1.93e-13         |" + LS
    + "                                 |G                0         |G                1         |G                2         |G                3         |G                4         |G                5         |G                6         |A                7         |T                8         |A                9         |A               10         |A               11         |A               12         |A               13         |A               14         |G               15         |G               16         |C               17         |G               18         |A               19         |C               20         |A               21         |T               22         |G               23         |C               24         |C               25         |A               26         |A               27         |T               28         |G               29         |T               30         |G               31         |T               32         |C               33         |G               34         |C               35         |C               36         |T               37         |T               38         |T               39         |T               40         |T               41         |C               42         |A               43         |A               44         |C               45         |T               46         |T               47         |T               48         |C               49         |C               50         |G               51         |A               52         |T               53         |" + LS
    ;
  public void testCGSpreadsheet() {
    final CgGotohEditDistance ed = new CgGotohEditDistance(7, new ScoreMatrixCGTest.MockRealignParamsCG(), 0);
    final String read = "        ATAAA AAGGCGACAT GCCAATGTGT        TTCAACTTTC";
    final String tmpl = "gggggggataAAA   AAGGCGACAT GCCAATGTGT CGCCTTTTTCAACTTTCCGATTAA";
    final byte[] s1 = DnaUtils.encodeString(read.replaceAll(" ", "").toLowerCase(Locale.ROOT));
    final byte[] s2 = DnaUtils.encodeString(tmpl.replaceAll(" ", "").toLowerCase(Locale.ROOT));
    //    final long start = System.nanoTime();
    final int[] v = ed.calculateEditDistance(s1, s1.length, s2, 7, Integer.MAX_VALUE, 7, true);
    //    final long end = System.nanoTime();
    //    System.out.println("time taken in seconds = " + (end - start) / 1.0E9);
    //    System.out.println(ed.toString());
    //    System.out.println("start: " + ActionsHelper.zeroBasedTemplateStart(v) + " score: " + ActionsHelper.alignmentScore(v));
    //    System.out.println("actions: " + ActionsHelper.toString(v));
    final AlignmentResult alignment = new AlignmentResult(s1, v, s2);
    assertEquals("1=1X21=7N10=", alignment.getCigarString(false, false));
    assertEquals(1, ActionsHelper.alignmentScore(v));
    assertEquals(10, ActionsHelper.zeroBasedTemplateStart(v));
    assertEquals(10, alignment.getStart());
    assertEquals(CG_SPREADSHEET, ed.toString());
  }

  public void testCGSpreadsheetRightArm() {
    final CgGotohEditDistance ed = new CgGotohEditDistance(7, new ScoreMatrixCGTest.MockRealignParamsCG(), 0);
    //final String read = "          ATAAA AAGGCGACAT GCCAATGTGT .....  TTCAACTTTC";
    //final String tmpl = "gggggggataAAA   AAGGCGACAT GCCAATGTGT CGCCTTTTTCAACTTTCCGATTAA";
    final String read = "       CTTTCAACTT         TGTGTAACCG TACAGCGGAA AAATA";
    final String tmpl = "AATTAGCCTTTCAACTT TTTCCGC TGTGTAACCG TACAGCGGAA   AAAataggggggg";
    final byte[] s1 = DnaUtils.encodeString(read.replaceAll(" ", "").toLowerCase(Locale.ROOT));
    final byte[] s2 = DnaUtils.encodeString(tmpl.replaceAll(" ", "").toLowerCase(Locale.ROOT));
    final int[] v = ed.calculateEditDistance(s1, s1.length, s2, 11, Integer.MAX_VALUE, 7, false);
    //    System.out.println(ed.toString());
    //    System.out.println("start: " + ActionsHelper.zeroBasedTemplateStart(v) + " score: " + ActionsHelper.alignmentScore(v));
    //    System.out.println("actions: " + ActionsHelper.toString(v));
    final AlignmentResult alignment = new AlignmentResult(s1, v, s2);
    assertEquals("10=7N21=1X1=", alignment.getCigarString(false, false));
    assertEquals(1, ActionsHelper.alignmentScore(v));
    assertEquals(7, ActionsHelper.zeroBasedTemplateStart(v));
    assertEquals(7, alignment.getStart());
    final String[] expectLines = CG_SPREADSHEET.split(LS);
    final String[] actualLines = ed.toString().split(LS);
    assertEquals(expectLines.length, actualLines.length);
    // should be all the same, except for the template numbering in the top and bottom lines.
    for (int i = 1; i < actualLines.length - 1; i++) {
      assertEquals(expectLines[i], expectLines[i], actualLines[i]);
    }
  }

  private void logContains(final String gap, final String range, final int... values) {
    final StringBuilder sb = new StringBuilder();
    sb.append(gap);
    sb.append(":\t");
    sb.append(range);
    for (final int val : values) {
      sb.append("\t");
      sb.append(val);
    }
    final String expected = sb.toString();
    assertTrue("Missing from log: " + expected + "\nLog: " + mLog, mLog.contains(expected));
  }

  public void testCGSnpInOverlap1() {
    checkCG("        ATAAG ACGGCGACAT GCCAATGTGT        TTCAACTTTC",
            "gggggggataAAA   ACGGCGACAT GCCAATGTGT CGCCTTTTTCAACTTTCCGATTAA",
            "ATAAGGGCGACATGCCAATGTGT.......TTCAACTTTC",
            "AAAACGGCGACATGCCAATGTGTCGCCTTTTTCAACTTTC",
            "=X==XBB====================NNNNNNN==========",
            "1=1X2=1X18=7N10=", 2, true, false, 7, true);
    //System.out.println(mLog);
    assertTrue(mLog.contains("CgGotohEditDistance mismatch=0.036180 delopen=0.001000 delextend=0.181694 insopen=0.002000 insextend=0.122598"));
    logContains("LeftOverlap", "-4..0", 0, 0, 1, 0, 0);
    logContains("LeftSmall", "0..3", 1, 0, 0, 0);
    logContains("LeftLarge", "4..8", 0, 0, 0, 1, 0);
    logContains("RightOverlap", "-4..0", 0, 0, 0, 0, 0);
    logContains("RightSmall", "0..3", 0, 0, 0, 0);
    logContains("RightLarge", "4..8", 0, 0, 0, 0, 0);
  }

  public void testCGSingleSnp() {
    checkCG("        AAAAC ACGGCGACAT GCCAtTGTGT        TTCAACTTTC",
        "gggggggattAAAAC   GGCGACAT GCCAATGTGT CGCCTT TTCAACTTTCCGATTAA",
        "AAAACGGCGACATGCCATTGTGT......TTCAACTTTC",
        "AAAACGGCGACATGCCAATGTGTCGCCTTTTCAACTTTC",
        "=====BB==============X=====NNNNNN==========",
        "17=1X5=6N10=", 1, true, false, 7, true);
    assertEquals(11, mSamRecord.getStartPosition());
  }
  public void testCGSingleSnpShift() {
    checkCG("        AAAAC ACGGCGACAT GCCAtTGTGT        TTCAACTTTC",
        "gggggggattAAAAC   GGCGACAT GCCAATGTGT CGCCTT TTCAACTTTCCGATTAA",
        "AAAACGGCGACATGCCATTGTGT......TTCAACTTTC",
        "AAAACGGCGACATGCCAATGTGTCGCCTTTTCAACTTTC",
        "=====BB==============X=====NNNNNN==========",
        "17=1X5=6N10=", 1, true, false, 5, true);
    assertEquals(11, mSamRecord.getStartPosition());
  }

  public void testCGLowestScore() {
    checkCG("        AAAAC ACGGCGACAT GCCAATGTGT        TTCAACTTTC",
        "gggggggattAAAAC   GGCGACAT GCCAATGTGT CGCCTT TTCAACTTTCCGATTAA",
        "AAAACGGCGACATGCCAATGTGT......TTCAACTTTC",
        "AAAACGGCGACATGCCAATGTGTCGCCTTTTCAACTTTC",
        "=====BB====================NNNNNN==========",
        "23=6N10=", 0, true, false, 7, true);
    assertEquals(11, mSamRecord.getStartPosition());
  }

  public void testCGSnpInOverlap2() {
    checkCG("        ATAAC AGGGCGACAT GCCAATGTGT        TTCAACTTTC",
          "gggggggataAAAAC   GGCGACAT GCCAATGTGT CGCCTTTTTCAACTTTCCGATTAA",
          "ATAACGGCGACATGCCAATGTGT.......TTCAACTTTC",
          "AAAACGGCGACATGCCAATGTGTCGCCTTTTTCAACTTTC",
          "=X===BB=X==================NNNNNNN==========",
          "1=1X21=7N10=", 2, true, false, 11, true);
    assertEquals(10, mAlignment.getStart());
  }

  public void testCGDelete1Mid() {
    checkCG("        ATAAA AAGGCGACAT GCCA TGTGTC        TTCAACTTTC",
            "gggggggataAAA   AAGGCGACAT GCCAATGTGTC GCCTTT TTCAACTTTCCGATTAA",
            "ATAAAGGCGACATGCCA TGTGTC......TTCAACTTTC".replaceAll(" ", ""),
            "AAAAAGGCGACATGCCAATGTGTCGCCTTTTTCAACTTTC",
            "=X===BB==============D======NNNNNN==========",
//            "| ||||||||||||||| ||||||      ||||||||||",
            "1=1X15=1D6=6N10=", 3, true, false, 7, true);
    assertEquals(10, mAlignment.getStart());
  }

  public void testCGDelete4Mid() {
    checkCG("        ATAAA AAGGCGACAT GCCA    GTCGCC        TTCAACTTTC",
            "gggggggataAAA   AAGGCGACAT GCCAATGTGTCGCC TTTTTT TTCAACTTTCCGATTAA",
            "ATAAAGGCGACATGCCA    GTCGCC......TTCAACTTTC".replaceAll(" ", ""),
            "AAAAAGGCGACATGCCAATGTGTCGCCTTTTTTTTCAACTTTC",
            "=X===BB==============DDDD======NNNNNN==========",
            "1=1X15=4D6=6N10=", 6, true, false, 7, true);
  }

  public void testCGInsert1Mid() {
    checkCG("        ATAAA AAGGCGACAT CCAATGTGTC        TTCAACTTTC",
        "GGGGGGGATAAAA   AAGGCGACAT CCA TGTGTC GCCTTT TTCAACTTTCCGATTAA",
        "ATAAAGGCGACATCCAATGTGTC......TTCAACTTTC",
        "AAAAAGGCGACATCCA-TGTGTCGCCTTTTTCAACTTTC",
        "=X===BB=============I======NNNNNN==========",
        "1=1X14=1I6=6N10=", 3, true, false, 7, true);
  }

  public void testCGInsert4Mid() {
    checkCG("        ATAAA AAGGCGACAT GATGTCTCGC        TTCAACTTTC",
        "gggggggataAAA   AAGGCGACAT     TCTCGC TTTTTT TTCAACTTTCCGATTAA",
        "ATAAAGGCGACATGATGTCTCGC......TTCAACTTTC",
        "AAAAAGGCGACAT----TCTCGCTTTTTTTTCAACTTTC",
        "=X===BB==========IIII======NNNNNN==========",
        "1=1X11=4I6=6N10=", 6, true, false, 7, true);
  }

  // overlap of 0 (cg_test_errors priors do allow this!)
  public void testCGOverlap0() {
    checkCG("acgtcacgtcacgtcacgtcacgtc     acgtcacgtc",
            "acgtcacgtcacgtcacgtcacgtcgggggacgtcacgtc",
            "ACGTCACGTCACGTCACGTCACGTC.....ACGTCACGTC",
            "ACGTCACGTCACGTCACGTCACGTCGGGGGACGTCACGTC",
            "=========================NNNNN==========",
            "25=5N10=", 0, false, true);
    assertEquals(0, mAlignment.getStart());
  }
  public void testCGOverlap0Rev() {
    checkCG("acgtcacgtcacgtcacgtcacgtc     acgtcacgtc",
            "gacgtgacgtcccccgacgtgacgtgacgtgacgtgacgt",
            "ACGTCACGTCACGTCACGTCACGTC.....ACGTCACGTC",
            "ACGTCACGTCACGTCACGTCACGTCGGGGGACGTCACGTC",
            "=========================NNNNN==========",
            "10=5N25=", 0, true, true);
    assertEquals(0, mAlignment.getStart());
  }

  // overlap of 4 (cg_test_errors priors do allow this!)
  public void testCGOverlap4() {
    checkCG("acgtc cgtc aacgtcacgtcacgtc     acgtcacgtc",
            "acgtc      aacgtcacgtcacgtcgggggacgtcacgtc",
            "ACGTCAACGTCACGTCACGTC.....ACGTCACGTC",
            "ACGTCAACGTCACGTCACGTCGGGGGACGTCACGTC",
            "=====BBBB====================NNNNN==========",
            "21=5N10=", 0, false, true);
    assertEquals(1, mSamRecord.getStartPosition());
    assertEquals(0, mAlignment.getStart());
  }

  // large gap of 8, small gap of 3, and overlap of 4 (right arm)
  public void testCG834RightArm() {
    checkCG("aacgtcacgt        cacgtcacgt   cacgtc cgtc cgtca",
            "aacgtcacgtggggggggcacgtcacgtgggcacgtc      cgtca",
            "AACGTCACGT........CACGTCACGT...CACGTCCGTCA",   //TODO: fix this: first 10 are "ccccccccccc"
            "AACGTCACGTGGGGGGGGCACGTCACGTGGGCACGTCCGTCA",
            "==========NNNNNNNN==========NNN==========BBBB=====",
            "10=8N10=3N11=", 0, false, false, 0, false);
    logContains("LeftOverlap", "-4..0", 0, 0, 0, 0, 0);
    logContains("LeftSmall", "0..3", 0, 0, 0, 0);
    logContains("LeftLarge", "4..8", 0, 0, 0, 0, 0);
    logContains("RightOverlap", "-4..0", 1, 0, 0, 0, 0);
    logContains("RightSmall", "0..3", 0, 0, 0, 1);
    logContains("RightLarge", "4..8", 0, 0, 0, 0, 1);
    assertEquals(1, mSamRecord.getStartPosition());
  }

  // large gap of 4, small gap of 0, and overlap of 0 (right arm)
  public void testCG400RightArm() {
    checkCG("aacgtcacgt     cacgtcacgtcacgtc cgtc cgtca",
            "aacgtcacgtgggg cacgtcacgtcacgtc cgtc cgtca",
            "AACGTCACGT....CACGTCACGTCACGTCCGTCCGTCA",
            "AACGTCACGTGGGGCACGTCACGTCACGTCCGTCCGTCA",
            "==========NNNN=========================",
            "10=4N25=", 0, false, false, 0, false);
    logContains("LeftOverlap", "-4..0", 0, 0, 0, 0, 0);
    logContains("LeftSmall", "0..3", 0, 0, 0, 0);
    logContains("LeftLarge", "4..8", 0, 0, 0, 0, 0);
    logContains("RightOverlap", "-4..0", 0, 0, 0, 0, 1);
    logContains("RightSmall", "0..3", 1, 0, 0, 0);
    logContains("RightLarge", "4..8", 1, 0, 0, 0, 0);
  }

  // overlap of 5 (cg_test_errors priors allow max of 4, so it chooses overlap of 4 plus an insert)
  // TODO: correct the cigar generation here (tricky because the overlap region contains a delete!)
  public void testCGOverlap5() {
    checkCG("      acgtc acgtc aacgtcacgtcacgt     acgtcacgtc",
            "ttgca acgtc       aacgtcacgtcacgtgggggacgtcacgtc",
            "ACGTCCAACGTCACGTCACGT.....ACGTCACGTC",
            "ACGT-CAACGTCACGTCACGTGGGGGACGTCACGTC",
            "====IBBBB====================NNNNN==========",
            "4=1I16=5N10=", 2, false, true);
  }

  public void testCGSmallGap2() {
    checkCG("acact ctgggggaaa    aacccccttt       acgtcacgtc",
            "acact   gggggaaa tt aacccccttt ggggg acgtcacgtc",
            "ACACTGGGGGAAA..AACCCCCTTT.....ACGTCACGTC",
            "ACACTGGGGGAAATTAACCCCCTTTGGGGGACGTCACGTC",
            "=====BB==========NN==========NNNNN==========",
            "13=2N10=5N10=", 0, false, true);
  }

  public void testCGSmallGap3() {
    checkCG("acact ctgggggaaa     aacccccttt       acgtcacgtc",
            "acact   gggggaaa ttt aacccccttt ggggg acgtcacgtc",
            "ACACTGGGGGAAA...AACCCCCTTT.....ACGTCACGTC",
            "ACACTGGGGGAAATTTAACCCCCTTTGGGGGACGTCACGTC",
            "=====BB==========NNN==========NNNNN==========",
            "13=3N10=5N10=", 0, false, true);
  }

  public void testCGSmallGap4() {
    // small gap of 4 is too big, so we get a gap of 3 plus 1 delete.
    checkCG("acact ctgggggaaa      aacccccttt       acgtcacgtc",
            "acact   gggggaaa tttt aacccccttt ggggg acgtcacgtc",
            "ACACTGGGGGAAA-...AACCCCCTTT.....ACGTCACGTC".replaceAll("-", ""),
            "ACACTGGGGGAAATTTTAACCCCCTTTGGGGGACGTCACGTC",
            "=====BB==========DNNN==========NNNNN==========",
            "13=1D3N10=5N10=", 2, false, true);
    logContains("LeftOverlap", "-4..0", 0, 0, 1, 0, 0);
    logContains("LeftSmall", "0..3", 0, 0, 0, 1);
    logContains("LeftLarge", "4..8", 0, 1, 0, 0, 0);
    logContains("RightOverlap", "-4..0", 0, 0, 0, 0, 0);
    logContains("RightSmall", "0..3", 0, 0, 0, 0);
    logContains("RightLarge", "4..8", 0, 0, 0, 0, 0);
  }

  public void testCGLargeGap4() {
    checkCG("acact ctgggggaaa aacccccttt      acgtcacgtc",
            "acact   gggggaaa aacccccttt gggg acgtcacgtc",
            "ACACTGGGGGAAAAACCCCCTTT....ACGTCACGTC",
            "ACACTGGGGGAAAAACCCCCTTTGGGGACGTCACGTC",
            "=====BB====================NNNN==========",
            "23=4N10=", 0, false, true);
    assertEquals(1, mSamRecord.getStartPosition());
  }

  public void testCGLargeGap7() {
    checkCG("acact ctgggggaaa aacccccttt         acgtcacgtc",
            "acact   gggggaaa aacccccttt ggggggg acgtcacgtc",
            "ACACTGGGGGAAAAACCCCCTTT.......ACGTCACGTC",
            "ACACTGGGGGAAAAACCCCCTTTGGGGGGGACGTCACGTC",
            "=====BB====================NNNNNNN==========",
            "23=7N10=", 0, false, true);
    assertEquals(1, mSamRecord.getStartPosition());
  }

  public void testCGLargeGap8() {
    checkCG("acact cactgggaaa    aacccccttt         acgtcacgtc",
            "acact     gggaaa ttt aacccccttt gggggggg acgtcacgtc",
            "ACACTGGGAAA...AACCCCCTTT........ACGTCACGTC",
            "ACACTGGGAAATTTAACCCCCTTTGGGGGGGGACGTCACGTC",
            "=====BBBB==========NNN==========NNNNNNNN==========",
            "11=3N10=8N10=", 0, false, true);
    logContains("LeftOverlap", "-4..0", 1, 0, 0, 0, 0);
    logContains("LeftSmall", "0..3", 0, 0, 0, 1);
    logContains("LeftLarge", "4..8", 0, 0, 0, 0, 1);
    logContains("RightOverlap", "-4..0", 0, 0, 0, 0, 0);
    logContains("RightSmall", "0..3", 0, 0, 0, 0);
    logContains("RightLarge", "4..8", 0, 0, 0, 0, 0);
  }

  public void testCGLargeGap9() {
    // large gap of 9 is too big, so we get a gap of 8 plus 1 delete.
    checkCG("acact ctgggggaaa aacccccttt           acgtcacgtc",
            "acact   gggggaaa aacccccttt ggggggggg acgtcacgtc",
            "ACACTGGGGGAAAAACCCCCTTT-........ACGTCACGTC".replaceAll("-", ""),
            "ACACTGGGGGAAAAACCCCCTTTGGGGGGGGGACGTCACGTC",
            "=====BB====================DNNNNNNNN==========",
            "23=1D8N10=", 2, false, true);
    assertEquals(1, mSamRecord.getStartPosition());
  }

  public void testCGComplex1() {
    // several MNPs and indels
    checkCG("acact ctgggACaaa  aac cctttgg        acgtcacgtc",
            "acact   gggggaaa aaacccctttgg gggggg acgtc cgtc",
            "ACACTGGGACAAA..AACCCTTTGG......ACGTCACGTC",
            "ACACTGGGGGAAAAAACCCCTTTGGGGGGGGACGTC-CGTC",
            "=====BB=====XX===NN=X========NNNNNN=====I====",
            "8=2X3=2N1=1X8=6N5=1I4=", 5, false, true);
  }

  public void testCGSoftClippingStart() {
    checkCG("tcaaa aaccggtttt gggggacgtc         cacgtcaagc",
            "  aaa   ccggtttt gggggacgtc actgctt cacgtcaagc",
            "TCAAACCGGTTTTGGGGGACGTC.......CACGTCAAGC",
            "NNAAACCGGTTTTGGGGGACGTCACTGCTTCACGTCAAGC",
            "XX===BB====================NNNNNNN==========",
            "2S21=7N10=", 2, false, true);
    assertEquals(-2, mAlignment.getStart());
  }

  public void testCGSoftClippingEnd() {
    checkCG("tcaaa aaccggtttt gggggacgtc         cacgtcaagc",
            "tcaaa   ccggtttt gggggacgtc actgctt cacgtca",
            "TCAAACCGGTTTTGGGGGACGTC.......CACGTCAAGC",
            "TCAAACCGGTTTTGGGGGACGTCACTGCTTCACGTCANNN",
            "=====BB====================NNNNNNN=======XXX",
            "23=7N7=3S", 3, false, true);
    assertEquals(1, mSamRecord.getStartPosition());
  }

  public void testCGSoftClippingBoth() {
    checkCG("tcaaa aaccggtttt gggggacgtc         cacgtcaagc",
            "  aaa   ccggtttt gggggacgtc actgctt cacgtca",
            "TCAAACCGGTTTTGGGGGACGTC.......CACGTCAAGC",
            "NNAAACCGGTTTTGGGGGACGTCACTGCTTCACGTCANNN",
            "XX===BB====================NNNNNNN=======XXX",
            "2S21=7N7=3S", 5, false, true);
    assertEquals(-2, mAlignment.getStart());
  }

  public void testCgInsertBeforeSpace() {
    checkCG("tatgagaaattctgggttgaaaatt     tttaagaatg",
          "gatatga  aattctgggttgaaaattcttttaagaatgttgaatattggcccccacggtcttctggcttgtagggtttctgcagag",
            "TATGAAATTCTGGGTTGAAAATT....TTTAAGAATG",
            "TATGAAATTCTGGGTTGAAAA--TTCTTTTAAGAATG",
            "=====BB==================IINNNN==========",
            "21=2I4N10=", 3, false, true);
    assertEquals(3, mSamRecord.getStartPosition());
  }

  public void testCgOverlapRealWorld6() {
    checkCG("attcttactanccccaactgagccc     agtaggagta",
            "atactcctacttttgctgggctcagttgggggaagtagaat",
            "ATTCTACTANCCCCAACTGAGCCC......AGTAGGAGTA",
            "ATTCTACTTCCCCCAACTGAGCCCAGCAAAAGTAGGAGTA",
            "=====B====XR==============NNNNNN==========",
            "10=6N14=2X8=", 1, true, true, 0, true, 0);
    assertEquals(2, mSamRecord.getStartPosition());
  }

  public void testCgOverlapRealWorld6NsMismatches() {
    checkCG("attcttactanccccaactgagccc     agtaggagta",
            "atactcctacttttgctgggctcagttgggggaagtagaat",
            "ATTCTACTANCCCCAACTGAGCCC......AGTAGGAGTA",
            "ATTCTACTTCCCCCAACTGAGCCCAGCAAAAGTAGGAGTA",
            "=====B====XR==============NNNNNN==========",
            "10=6N14=2X8=", 2, true, true, 0, true, 1);
    assertEquals(2, mSamRecord.getStartPosition());
  }

  public void testNgsCgExample() {
    checkCG("TGCATGCATGCATGCACAAAGCTAC      ATTGCTATGC",
            "TGCATGCATGCATGCACATTGCTACAAAAAAATTGCTATGCCCATTACG",
            "TGCATGCATGCATGCACAAAGCTAC......ATTGCTATGC",
            "TGCATGCATGCATGCACATTGCTACAAAAAAATTGCTATGC",
            "==================XX=====NNNNNN==========",
            "18=2X5=6N10=", 2, true, false, 2, true);
  }

  public void testValidation() {
    checkCG("GTACTGTGTGCGTCCTTGACATCTA     ACCGGCGCCT",
            "GTAGT  GTGCGTCCTTGACATCTAAGGATATCGGCGCCTGA",
            "GTACTGTGCGTCCTTGACATCTA.....ACCGGCGCCT",
            "GTAGTGTGCGTCCTTGACATCTAAGGATATCGGCGCCT",
            "===X=BB====================NNNNN=X========",
            "3=1X19=5N1=1X8=", 2, true, false, 0, true);
  }


  /* This test is dependent on the cg errors priors :/
   public void testBLah() {
    checkCG("TCCATTAAAC.....TTCTTTATAAATTATCCTTCTTTTC",
            "ATGAGTCCATTAAACCTTTTATTCTTTATAAATTATCCAGTCACAGATATTTCTTC",
                 "tccattaaac......ttctttataaattatcct-ttttc",
                 "tccattaaaccttttattctttataaattatccag-tcac",
                 "||||||||||      |||||||||||||||||   |  |",
            "10=6N17=1X1X1=2X1=", 5, false, false, 5, false);

  }*/

  public void testOverlapWeirdScoring() {
    checkCG("tcaaa aagcggtttt gggggacgtc         cacgtcaagc",
            "tcatt ccggtttt gggggacgtc actgctt cacgtcaagcatatat",
            "TCAAAGCGGTTTTGGGGGACGTC.......CACGTCAAGC",
            "TCATTCCGGTTTTGGGGGACGTCACTGCTTCACGTCAAGC",
            "===XXBBXXX=================NNNNNNN==========",
            "3=3X17=7N10=", 5, false, true);
    assertEquals(1, mSamRecord.getStartPosition());    //remember, -4 is 1-based
  }

  public void testLargeSoftClipLeft() {
    checkCG("tcaaa aaccggtttt gggggacgtc         cacgtcaagc",
            "        ccggtttt gggggacgtc actgctt cacgtcaagcatatat",
            "TCAAACCGGTTTTGGGGGACGTC.......CACGTCAAGC",
            "NNNNNCCGGTTTTGGGGGACGTCACTGCTTCACGTCAAGC",
            "XXXXXBBXX==================NNNNNNN==========",
            "5S18=7N10=", 7, false, true);
    assertEquals(-5, mAlignment.getStart());
  }
  public void testLargeSoftClipRight() {
    checkCG("aaaaa aaccggtttt gggggacgtc         cacgtcaagc",
            "aaaaa aaccggtttt gggggacgtc actgctt cacgt",
            "AAAAACCGGTTTTGGGGGACGTC.......CACGTCAAGC",
            "AAAAACCGGTTTTGGGGGACGTCACTGCTTCACGTNNNNN",
            "=====BB====================NNNNNNN=====XXXXX",
        "23=7N5=5S", 5, false, true);
    assertEquals(3, mSamRecord.getStartPosition());    //remember, -4 is 1-based
  }
  public void testLargeSoftClipRight2() {
    checkCG("      cgaactgcac         ctgcaggggg ttttggccaa aaact",
            "tatatacgaactgcac ttcgtca ctgcaggggg ttttggcc",
            "CGAACTGCAC.......CTGCAGGGGGTTTTGGCCAAACT",
            "CGAACTGCACTTCGTCACTGCAGGGGGTTTTGGCCNNNNN",
            "==========NNNNNNN==================XXBBXXXXX",
            "10=7N18=5S", 7, false, false, 7, false);
    assertEquals(7, mSamRecord.getStartPosition());
  }

  public void testEvenLargerSoftLeftClip() {
    checkCG("tcaaa aaccggtttt gggggacgtc         cacgtcaagc",
            "          ggtttt gggggacgtc actgctt cacgtcaagcatatat",
            "TCAAACCGGTTTTGGGGGACGTC.......CACGTCAAGC",
            "NNNNNNNGGTTTTGGGGGACGTCACTGCTTCACGTCAAGC",
            "XXXXXBBXXXX================NNNNNNN==========",
            "7S16=7N10=", 9, true, false, -8, true);
    assertEquals(-7, mAlignment.getStart());
  }

  public void testCgOverlapRealWorld7() {
      checkCG("tcaanaacanggaacagagcagggg     acttgcagaa",
              "ttctgcaagtgttgggcccctgctctgttccctgtttgatatactggggg",
              "TCAANCANGGAACAGAGCAGGGG......ACTTGCAGAA",
              "TCAAACAGGGAACAGAGCAGGGGCCCAACACTTGCAGAA",
              "====RBB====R===============NNNNNN==========",
              "10=6N15=1X2=1X4=", 0, true, true);
  }

  static final String TEMPLATE = ""
    + "GGTGGTCGAGTTGTGTGTAGGCGATGGCACGACCTCAGGAGTTACTGTACCATGGTTGCTTCTAGCGAGAAATTCTTG"
    ;
  static final String[] READS = {
    //One-Based-Start,  Read-String,                  Cigar
    "6",  "TCGAGAGTTGTGTGTAGGCGATGGC     TCAGGAGTTA".replaceAll(" ", ""), "23=6N10=",
    "16", "TGTAGAGGCGATGGCACCTCAGGAG     ACCATGGTTG".replaceAll(" ", ""), "13=3N10=7N10=",
    "21", "GCGATGATGGCACGACCTCAAGAGT     ACCATGGTTG".replaceAll(" ", ""), "17=1X4=6N10=",
    "22", "CGATGTGGCACGACCTCAGGAGTTA     CATGGTTGCT".replaceAll(" ", ""), "23=6N10=",
    "40", "AGTTATTCTGTACCATGGTTGCTTC     AATTCTTGCA".replaceAll(" ", ""), "23=8N8=2S"
  };
  public void testMultiReads() throws InvalidParamsException, IOException {
    final ByteArrayOutputStream bos = new ByteArrayOutputStream();
    Diagnostic.setLogStream(new PrintStream(bos));
    final RealignParamsImplementation params = new RealignParamsImplementation(MachineErrorParams.builder().errors("cg_test_errors").create());
    final CgGotohEditDistance ed = new CgGotohEditDistance(7, params, 0);
    final byte[] tmpl = DnaUtils.encodeString(TEMPLATE);
    for (int i = 0; i < READS.length; i += 3) {
      final int start = Integer.parseInt(READS[i]);
      final byte[] read = DnaUtils.encodeString(READS[i + 1]);
      final String cigar = READS[i + 2];
      final int approxStart = start ^ 5; // move the start around a bit
      final int[] v = ed.calculateEditDistance(read, read.length, tmpl, approxStart, Integer.MAX_VALUE, 7, true);
      //      System.out.println(ed.toString());
      //      System.out.println("start: " + ActionsHelper.zeroBasedTemplateStart(v)
      //          + " score: " + ActionsHelper.alignmentScore(v));
      //      System.out.println("actions: " + ActionsHelper.toString(v));
      final AlignmentResult alignment = new AlignmentResult(read, v, tmpl);
      alignment.setIdentifyingInfo(true, false);
      assertEquals(cigar, alignment.getCigarString(false, false));
      assertEquals(start - 1, ActionsHelper.zeroBasedTemplateStart(v));
      assertEquals(start - 1, alignment.getStart());
    }
    ed.logStats();
    mLog = bos.toString();
    logContains("LeftOverlap", "-4..0", 0, 1, 4, 0, 0);
    logContains("LeftSmall", "0..3", 4, 0, 0, 1);
    logContains("LeftLarge", "4..8", 0, 0, 3, 1, 1);
    logContains("RightOverlap", "-4..0", 0, 0, 0, 0, 0);
    logContains("RightSmall", "0..3", 0, 0, 0, 0);
    logContains("RightLarge", "4..8", 0, 0, 0, 0, 0);
    Diagnostic.setLogStream();
  }

  public static void main(String[] args) {
    new CgGotohEditDistanceTest().testCGLeftRightSymmetry();
  }

  public void testCGLeftRightSymmetry() {
    getEditDistanceInstance(1, 1, 0); // new CgGotohEditDistance(7, new ScoreMatrixCGTest.MockRealignParamsCG());
    final CgGotohEditDistance ed = mCgGotoh;
    final String[] reads = {
        "ATAAA AAGGCGACAT GCCAATGTGT        TTCAACTTTC",
        "ACGTT CGTTACGTAC GTACGTACGT        ACGACGACGA",
        "GGGAA ACATATATAT GGCGACATGC        ATTTTTACAC",
        "GGGTG TGCAATTGAC ACGACGACGA        TTCAACTTTC",
        "ATCCA ACGGCCACAT GCCATTTTGT        CGTTACGTAC",
        "AAAAA AAAAAAAAAA AAAAAAAAAA        AAAAAAAAAA",
        "AAAAA AACACACACA CACACACACA        ACACACACAA",
        "AAAAA CCCCCCCCCC GGGGGGGGGG        TTTTTTTTTT",
        "GACGT ACGTATACGT GGTGACTTGC        ATACCCACGC",
        "GACGT GCATACTATC GTATAAAAGC        TCGCCCACTT",
    };
    final PortableRandom rand = new PortableRandom(42);
    final String[] mutations = {"A", "C", "G", "T", "", "AC", "GT"};
    String mut = "";
    final int[] hist = new int[40];
    for (int repeat = 0; repeat < reads.length; repeat++) {
      final String read = reads[repeat % reads.length].replaceAll(" ", "").toLowerCase(Locale.ROOT);
      String tmpl = "aaaagggg" + read.substring(0, 25) + "tatata" + read.substring(25) + "accaacca";
      //System.out.println("@@@@@@@@@@@@@@@ starting with read:" + read + "@@@@@@@@@@@@@@@@@@@");
      int readStart = 8;
      int readEnd = readStart + 39;  // ie. 35 - 2 + 0 + 6 (for the most common gap sizes)
      for (int age = 0; age < 10 && Math.abs(readEnd - readStart - 39) < 5 ; age++) {
        final byte[] s1 = DnaUtils.encodeString(read);
        final byte[] s2 = DnaUtils.encodeString(tmpl.toLowerCase(Locale.ROOT));
        final int[] tmp = ed.calculateEditDistance(s1, s1.length, s2, readStart, Integer.MAX_VALUE, 7, true);
        final int[] v = new int[tmp.length];
        System.arraycopy(tmp, 0, v, 0, tmp.length);
        //      System.out.println(ed.toString());
        //      System.out.println("start: " + ActionsHelper.zeroBasedTemplateStart(v)
        //          + " score: " + ActionsHelper.alignmentScore(v));
        //      System.out.println("actions: " + ActionsHelper.toString(v));
        final int leftScore = ActionsHelper.alignmentScore(v);

        // now check that we get the same alignment score if we flip read+template onto the opposite arm.
        final byte[] s1R = DnaUtils.encodeString(reverse(read).replaceAll(" ", "").toLowerCase(Locale.ROOT));
        final byte[] s2R = DnaUtils.encodeString(reverse(tmpl).replaceAll(" ", "").toLowerCase(Locale.ROOT));
        assertEquals(s1.length, s1R.length);
        assertEquals(s2.length, s2R.length);
        final int startR = s2.length - readEnd ;
        final int[] vR = ed.calculateEditDistance(s1R, s1R.length, s2R, startR, Integer.MAX_VALUE, 7, false);
        final int rightScore = ActionsHelper.alignmentScore(vR);
        //System.out.println(mCgGotoh.toString());
        //assertTrue(Math.abs(leftScore - rightScore) <= 1);  // but they should be equal!
        final int diff = leftScore - rightScore;
        final String actsLeft = ActionsHelper.toString(v);
        final String actsRight = ActionsHelper.toString(vR);
        final StringBuilder sb = new StringBuilder(actsRight);
        sb.reverse();
        hist[diff + 20]++;
        if (diff != 0) {
          System.out.println("read=" + StringUtils.getSpaceString(readStart + 1) + read);
          System.out.println("tmpl=" + tmpl.substring(0, readStart) + " " + tmpl.substring(readStart, readEnd) + " " + tmpl.substring(readEnd));
          System.out.println("" + actsLeft + " " + leftScore + " " + rightScore + " " + sb.toString() + "  last mutation: " + mut + " age=" + age);
        }
        assertEquals(leftScore, rightScore);
        assertEquals(actsLeft, sb.toString());
        // now mutate the template a bit
        final int r = rand.nextInt(mutations.length * (tmpl.length() - 2));
        mut = mutations[r % mutations.length];
        final int pos = r / mutations.length;
        if (pos <= readStart) {
          readStart = readStart - 1 + mut.length();
        }
        if (pos <= readEnd) {
          readEnd = readEnd - 1 + mut.length();
        }
        tmpl = tmpl.substring(0, pos) + mut + tmpl.substring(pos + 1);
      }
    }
    //System.out.println("Score Differences: " + java.util.Arrays.toString(hist));
  }

  /* TODO: add in some of these tests too.  (They usually need some modification).
  public void testOverlap2Ugly() {

      checkCG("tatccgaaattctgggttgaaaatt     tttaagaatg",
              "gatatga  aattctgggttgaaaatt  cttttaagaatgttgaatattggcccccacggtcttctggcttgtagggtttctgcagag",
              "tatccgaaattctgggttgaaaatt.....tttaagaatg",
              "tat-.gaaattctgggttgaaaa--ttctt-ttaagaatg",
              "|||  ||||||||||||||||||        |||||||||",
              "3=1I18=2I5N1I9=", 7, false, true);
  }
  public void testOverlap2() {
      check("tatgagaaattctgggttgaaaatt     tttaagaatg",
            "gatatgaaattctgggttgaaaattcttttaagaatgttgaatattggcccccacggtcttctggcttgtagggtttctgcagag",
            "tatgagaaattctgggttgaaaatt.....tttaagaatg",
            "tat..gaaattctgggttgaaaa--ttctt-ttaagaatg",
            "|||  ||||||||||||||||||        |||||||||",
            "21=2I5N1I9=", 5, false);
  }
  public void testOverlap2RC() {
      check("tatgagaaattctgggttgaaaatt     tttaagaatg",
            "ccaatattcaacattcttaaaagaattttcaacccagaatttcatatc",
            "tatgagaaattctgggttgaaaatt.....tttaagaatg",
            "tat..gaaattctgggttgaaaa--ttctt-ttaagaatg",
            "|||  ||||||||||||||||||        |||||||||",
            "9=1I5N2I21=", 5, true);
  }

  public void testOverlap2Rev() {
      check("gtaagaattt     ttaaaagttgggtcttaaagagtat",
            "gtaagaattttcttaaaagttgggtcttaaagtatag",
            "gtaagaattt.....ttaaaagttgggtcttaaagagtat",
            "gtaagaa---ttttcttaaaagttgggtcttaaag..tat",
            "|||||||        ||||||||||||||||||||  |||",
            "7=3I5N23=", 4, false);
  }

  public void testCgOverlap1() {
      check("tatcaaaaattctgggttgaaaatt     tttaagaatg",
            "gatatcaaaattctgggttgaaaattcttttaagaatgttg",
            "tatcaaaaattctgggttgaaaatt.....tttaagaatg",
            "tatc.aaaattctgggttgaaaa--ttctt-ttaagaatg",
            "|||| ||||||||||||||||||        |||||||||",
            "22=2I5N1I9=", 5, false);
  }

  public void testCgOverlap1Rev() {
      check("gtaagaattt     ttaaaagttgggtcttaaaaactat",
            "gttgtaagaattttcttaaaagttgggtcttaaaactatag",
            "gtaagaattt.....ttaaaagttgggtcttaaaaactat",
            "gtaagaa---ttttcttaaaagttgggtcttaaaa.ctat",
            "|||||||        |||||||||||||||||||| ||||",
            "7=3I5N24=", 4, false);
  }
  public void testCgOverlap1RevRC() {
      check("gtaagaattt     ttaaaagttgggtcttaaaaactat",
            "ctatagttttaagacccaacttttaagaaaattcttacaac",
            "gtaagaattt.....ttaaaagttgggtcttaaaaactat",
            "gtaagaa---ttttcttaaaagttgggtcttaaaa.ctat",
            "|||||||        |||||||||||||||||||| ||||",
            "24=5N3I7=", 4, true);
  }

  //  not supporting overlap of 3 properly yet.
  public void testCgOverlap3() {
    final AlignmentResult ar =
      check("atgaagaaattctgggttgaaaatt     tttaagaatg",
            "gatatgaaattctgggttgaaaattcttttaagaatgttgaatat",
            "atgaagaaattctgggttgaaaatt.....tttaagaatg",
            "at---gaaattctgggttgaaaa--ttctt-ttaagaatg",
            "||   ||||||||||||||||||        |||||||||", 5, false, 3);
    final String cigar = CigarFormatter.alignmentToCigar(ar.mRead, ar.mTemplate, TOTAL_REFLENGTH_LONG, ar.mMatch, false, ar.mStart);
    assertEquals("2M3I18M2I5N1I9M", cigar);
  }

  // Test that '.' in a read is treated the same as 'n'.
  public void testAlignWithIndelsCG5() {
    final String read = "tccttag     gatgcaagacaagagggcctctc";
    final byte[] rr = DnaUtils.encodeArray(read.getBytes());
    EditDistance ed = getEditDistanceInstance(1, 1);
    final String t = "aaaaaatccttagactcagaggatgcaagacaagagggcctctcgaatgt";
                          //tccttag   .....gatgcaagacaagagggcctctc
    final byte[] template = DnaUtils.encodeArray(t.getBytes());
    int[] actions = ed.calculateEditDistance(rr, rr.length, template, 7, false, Integer.MAX_VALUE, 9);
    final AlignmentResult ar = new AlignmentResult(rr, actions, template.length, template, false);

    ar.setIdentifyingInfo(true, false);
    ar.setRemainingOutput(1, "toxic");
    final SAMRecord sr = ar.toSamRecord(header(), false, null);
    final String srformatted = sr.format();
    assertTrue(srformatted.contains("TCCTTAGGATGCAAGACAAGAGGGCCTCTC"));
    assertTrue(srformatted.contains("7=7N1D23="));
    assertTrue(srformatted.contains("AS:i:2\tNM:i:1"));

    final String rcTemplate = DnaUtils.reverseComplement(t);
    ed = getEditDistanceInstance(1, 1);
    final byte[] rctemplateBytes = DnaUtils.encodeArray(rcTemplate.getBytes());
    actions = ed.calculateEditDistance(rr, rr.length, rctemplateBytes, 7, true, Integer.MAX_VALUE, 9);
    final AlignmentResult ar2 = new AlignmentResult(rr, actions, rctemplateBytes.length, rctemplateBytes, false);

    ar2.setIdentifyingInfo(true, true);
    ar2.setRemainingOutput(1, "toxic");
    final SAMRecord sr2 = ar2.toSamRecord(header(), false, null);

    final String sr2formatted = sr2.format();
    assertTrue(sr2formatted.contains("GAGAGGCCCTCTTGTCTTGCATCCTAAGGA"));
    assertTrue(sr2formatted.contains("23=1D7N7="));
    assertTrue(sr2formatted.contains("AS:i:2\tNM:i:1"));
  }

  public void testCgOverlapRealWorld() {
      check("ataaaaaggcgacatgccaatgtgt     ttcaactttc",
            "ataaaaaaggcgacatgccaatgtgtcgcctttttcaactttccgattaagaacctgctcagcgggt",
            "ataaaaaggcgacatgccaatgtgt.......ttcaactttc",
            "aaa..aaggcgacatgccaatgtgtcgcctttttcaactttc",
            "| |  ||||||||||||||||||||       ||||||||||",
            "1=1X21=7N10=", 1, false);
  }
  public void testCgOverlapRealWorld2() {
      check("taacgccgat     ccaacggttatctcgatttttttta",
            "taacgccgattgaggccaacggttatctcgatttttttatc",
            "taacgccgat.....ccaacggttatctcgatttttttta",
            "taacgccgattgaggccaacggttatctcgatttt.ttta",
            "||||||||||     |||||||||||||||||||| ||||",
            "10=5N24=", 0, false);
  }

  public void testCgOverlapRealWorld3() {
      check("cgacggcggt     tgcttcactatgggaaagaggtcca",
            "cgacggcggtgggattgcttcactatgggaaagagtcacat",
            "cgacggcggt.....tgcttcactatgggaaagaggtc-ca",
            "cgacggcggtgggattgcttcactatgggaaagag.tcaca",
            "||||||||||     |||||||||||||||||||| || ||",
            "10=5N22=1D2=", 2, false);
  }

  public void testCgOverlapRealWorld4() {
      check("tttgtgtaggtcggataaggcgttc     atccgacacg",
            "ttttgtaggtcggataaggcgttcacgccgcatccgacacgg",
            "tttgtgtaggtcggataaggcgttc.......atccgacacg",
            "ttt..gtaggtcggataaggcgttcacgccgcatccgacacg",
            "|||  ||||||||||||||||||||       ||||||||||",
            "23=7N10=", 0, false, 1);
  }

  public void testCgOverlapRealWorld5() {
    final AlignmentResult ar =
      check("cgacggcggt     tgcttcactatgggaaagaggtcca",
            "cgacggcggtgggattgcttcactatgggaaagagtcacat",
            "cgacggcggt.....tgcttcactatgggaaagaggtc-ca",
            "cgacggcggtgggattgcttcactatgggaaagag.tcaca",
            "||||||||||     |||||||||||||||||||| || ||",
            "10=5N22=1D2=", 2, false, 1);
    assertEquals(1, ar.mismatches());
  }



  public void testCgRegression() {
    final AlignmentResult ar =
      check("cnggttaaaatatgaagtgaccacc     atgcttgaga",
            "ctggtaaaatatgaagtgaccaccaaagggagcttgagagaggagaaaatgact",
            "cnggttaaaatatgaagtgaccacc.....atgcttgaga",
            "ctgg.taaaatatgaagtgaccaccaaagggagcttgaga",
            "|||| ||||||||||||||||||||       ||||||||",
            "24=5N2X8=", 2, false, 1);
    assertEquals(2, ar.mismatches());
  }

  public void testFromAlignmentResult() {
    final AlignmentResult ar =
    check("tagacaaatg     ggattacaagccacaggagggggaa",
          "tagacaaatgtgactggattacaagccacaggaga",
          "tagacaaatg.....ggattacaagccacaggagggggaa",
          "tagacaaatgtgactggattacaagccacaggaga..nnn",
          "||||||||||     |||||||||||||||||||      ",
            "10=5N19=1X3S", 4, false, 1);
    assertEquals(4, ar.mismatches());
  }

  public void testFromAlignmentResult2() {
    final AlignmentResult ar =
    check("ctgctgaccc     gacctactggcggactcgggggtta",
            "ctgctgacccagagccggacctactggcggactcgggtta",
            "ctgctgaccc.......gacctactggcggactcgggggtta",
            "ctgctgacccagagccggacctactggcggactcggg..tta",
            "||||||||||       ||||||||||||||||||||  |||",
            "10=7N23=", 0, false, 1);
    assertEquals(0, ar.mismatches());
  }

  public void testBug609SoftClippingEnd() {
    final AlignmentResult ar =
      check("tagacaaatg     ggattacaagccacaggagtgaata",
            "tagacaaatgtgactggattacaagccacaggaga",
            "tagacaaatg.....ggattacaagccacaggagtgaata",
            "tagacaaatgtgactggattacaagccacaggaga.nnnn",
            "||||||||||     |||||||||||||||||||      ",
            "10=5N19=1X4S", 5, false);
    assertEquals(5, ar.mismatches());
  }

  public void testBug609SoftClippingStart() {
    final AlignmentResult ar =
      check("ttctgtggggtgactggattacaag     ggagtgaata",
            "tggggtgactggattacaagaaaaaggagtgaat",
            "ttctgtggggtgactggattacaag.....ggagtgaata",
            "nnn..tggggtgactggattacaagaaaaaggagtgaatn",
            "     ||||||||||||||||||||     ||||||||| ",
            "3S20=5N9=1X", 4, false, -3);
    assertEquals(4, ar.mismatches());
  }

  private AlignmentResult getAlignment(final String read, final String template, final EditDistance f, final boolean rc) {
    return getAlignment(read, template, f, 0, rc);
  }

  private AlignmentResult getAlignment(final String read, final String template, final EditDistance f, final int startPos, final boolean rc) {
    final byte[] s1 = read.getBytes();
    final byte[] s2 = template.getBytes();
    DnaUtils.encodeArray(s1);
    DnaUtils.encodeArray(s2);
    final int[] actions = f.calculateEditDistance(s1, s1.length, s2, startPos, rc, Integer.MAX_VALUE, 9);
    final AlignmentResult alignment = new AlignmentResult(s1, actions, s2.length, s2, false);
    if (rc) {
      alignment.setIdentifyingInfo(alignment.isFirst(), true);
    }
    return alignment;
  }
  public void testMultipleCgAlignments() {
    // Check that any object churn reduction isn't overriding something.
    final EditDistance f = getEditDistanceInstance(1, 1);
    final AlignmentResult alignment3 = getAlignment("attcttactanccccaactgagccc     agtaggagta",
            "atactcctacttttgctgggctcagttgggggaagtagaat", f, true);
    final String read1 = "cgacggcggt.....tgcttcactatgggaaagaggtcca";
    final AlignmentResult alignment1 = getAlignment(read1, "cgacggcggtgggattgcttcactatgggaaagagtcacat", f, false);
    final String template = "ctatagttttaagacccaacttttaagaaaattcttacaac";
    final String read2 = "gtaagaattt.....ttaaaagttgggtcttaaaaactat";
    final AlignmentResult alignment2 = getAlignment(read2, template, f, true);
    alignment2.setIdentifyingInfo(true, true);
    final String cigar1 = alignment1.getCigarString(alignment1.getStart(), false);
    final String cigar2 = alignment2.getCigarString(alignment2.getStart(), true);
    assertEquals("10=5N22=1D2=", cigar1);
    assertEquals("24=5N3I7=", cigar2);
    assertEquals("gtaagaattt.....ttaaaagttgggtcttaaaaactat" + "\t"
            + "gtaagaa---ttttcttaaaagttgggtcttaaaa.ctat" + "\t"
            + "|||||||        |||||||||||||||||||| ||||" , alignment2.tabularString());
    assertEquals("cgacggcggt.....tgcttcactatgggaaagaggtc-ca" + "\t"
            + "cgacggcggtgggattgcttcactatgggaaagag.tcaca" + "\t"
            + "||||||||||     |||||||||||||||||||| || ||", alignment1.tabularString());
    assertEquals("attcttactanccccaactgagccc......agtaggagta" + "\t"
            + "attc.tacttcccccaactgagcccagcaaaagtaggagta" + "\t"
            + "|||| |||| |||||||||||||||      ||||||||||", alignment3.tabularString());
  }
   */

  /* public void testRcBullcrap() {

    final EditDistance ed = new CgEditDistance(1);
    final AlignmentResult ar2 = getAlignment("attaaaaaaattttttttttttttt.....gacaaagtct", "cgacagagcgagattccgtctcaaacaaaaaaaaaaaaaaaaaaaattagctgggtgtgattataggtgc", ed, 16, true);
      //                        gcacctataatcacacccagctaattttttttttttttttttttgtttgagacggaatctcgctctgtcg < 14 start                           x
      //                   gcacctataatcacacccagctaattttttttttttttttttttgtttgagacggaatctcgctctgtcg < 19 start
    System.err.println(ar2.getScore());
    System.err.println(ar2.getActionsString());

    System.err.println("-------------------------------------------------");
    System.err.println("-------------------------------------------------");
    System.err.println("-------------------------------------------------");
    System.err.println("-------------------------------------------------");

    final EditDistance ed2 = new LoopingEditDistance(new RcEditDistance(new CgEditDistance(1)));
    final AlignmentResult ar3 = getAlignment("attaaaaaaattttttttttttttt.....gacaaagtct", "cgacagagcgagattccgtctcaaacaaaaaaaaaaaaaaaaaaaattagctgggtgtgattataggtgc", ed2, 16, true);
    System.err.println(ar3.getScore());
    System.err.println(ar3.getActionsString());
    assertEquals(8, ar3.getScore());

    System.err.println("-------------------------------------------------");
    System.err.println("-------------------------------------------------");
    System.err.println("-------------------------------------------------");
    System.err.println("-------------------------------------------------");


    final EditDistance ed3 = new RcEditDistance(new LoopingEditDistance(new CgEditDistance(1)));
    final AlignmentResult ar4 = getAlignment("attaaaaaaattttttttttttttt.....gacaaagtct", "cgacagagcgagattccgtctcaaacaaaaaaaaaaaaaaaaaaaattagctgggtgtgattataggtgc", ed3, 16, true);
    System.err.println(ar4.getScore());
    System.err.println(ar4.getActionsString());

    assertEquals(8, ar4.getScore());
  }*/

  public void testFormat() {
    assertEquals("         ", CgGotohEditDistance.format(-0.0));
    assertEquals("         ", CgGotohEditDistance.format(0.0));
    assertEquals(" 1.00e+00", CgGotohEditDistance.format(1.0));
    assertEquals(" 5.56e-01", CgGotohEditDistance.format(0.55555));
    assertEquals(" 9.99e-11", CgGotohEditDistance.format(0.0000000000999));
    assertEquals("     -inf", CgGotohEditDistance.format(Double.NEGATIVE_INFINITY));
    assertEquals("      inf", CgGotohEditDistance.format(Double.POSITIVE_INFINITY));
    assertEquals("      NaN", CgGotohEditDistance.format(Double.NaN));
  }

  public void testRowOffsetLeftArm() {
    getEditDistanceInstance(1, 1, 0);
    final int[] offsets = mCgGotoh.makeRowOffsets(CgGotohEditDistance.LEFT_ARM);
    assertEquals(36, offsets.length);
    final int half = -8; // half of the width of the band (at row 1).
    assertEquals(half, offsets[0]);
    assertEquals(half + 1, offsets[1]);
    assertEquals(half + 5, offsets[5]);
    assertEquals(half + 6 - 2, offsets[6]);
    assertEquals(half + 15 - 2, offsets[15]);
    assertEquals(half + 16 - 2, offsets[16]);
    assertEquals(half + 25 - 2, offsets[25]);
    assertEquals(half + 26 - 2 + 6, offsets[26]);
    assertEquals(half + 35 - 2 + 6 , offsets[35]);
  }

  public void testRowOffsetRightArm() {
    getEditDistanceInstance(1, 1, 0);
    final int[] offsets = mCgGotoh.makeRowOffsets(CgGotohEditDistance.RIGHT_ARM);
    assertEquals(36, offsets.length);
    final int half = -8; // half of the width of the band (at row 1).
    assertEquals(half, offsets[0]);
    assertEquals(half + 1, offsets[1]);
    assertEquals(half + 10, offsets[10]);
    assertEquals(half + 11 + 6, offsets[11]);
    assertEquals(half + 20 + 6, offsets[20]);
    assertEquals(half + 21 + 6, offsets[21]);
    assertEquals(half + 30 + 6, offsets[30]);
    assertEquals(half + 31 + 6 - 2, offsets[31]);
    assertEquals(half + 35 + 6 - 2, offsets[35]);
  }

  public void testFixed() {
    getEditDistanceInstance(1, 1, 0);
    assertNull(mCgGotoh.calculateEditDistanceFixedStart(null, 0, 0, null, 0, 0, 0));
    assertNull(mCgGotoh.calculateEditDistanceFixedEnd(null, 0, 0, null, 0, 0, 0, 0));
    assertNull(mCgGotoh.calculateEditDistanceFixedBoth(null, 0, 0, null, 0, 0, 0, 0));
  }

  public void testSoftClip() {
    final MemoryPrintStream mps = new MemoryPrintStream();
    Diagnostic.setLogStream(mps.printStream());
    try {
      getEditDistanceInstance(1, 1, 0);

      final byte[] t = DnaUtils.encodeString("gatcatcggcgacatgccattgtgttttttttcaacttt");
      final byte[] s1 = DnaUtils.encodeString("gatcatcggcgacatgccattgtgt     ttcaactttc".replaceAll(" ", ""));
      final byte[] t2 = DnaUtils.encodeString("gatcatcggcgacatgccattgtgttttttttcaa");

      //System.err.println(s1.length);
      int[] actions = mCgGotoh.calculateEditDistance(s1, s1.length, t, 0, 10, 5, true);
      if (actions != null) {
        assertEquals("=========================NNNNN=========X", ActionsHelper.toString(actions));
        assertEquals(1, actionsAlignmentScore(actions, false));
        assertEquals(0, actions[ActionsHelper.TEMPLATE_START_INDEX]);
      }

      actions = mCgGotoh.calculateEditDistance(s1, s1.length, t2, 0, 10, 5, true);
      if (actions != null) {
        assertEquals("=========================NNNNN=====XXXXX", ActionsHelper.toString(actions));
        assertEquals(5, actionsAlignmentScore(actions, false));
        assertEquals(0, actions[ActionsHelper.TEMPLATE_START_INDEX]);
      }

      getEditDistanceInstance(1, 1, 1);

      actions = mCgGotoh.calculateEditDistance(s1, s1.length, t, 0, 10, 5, true);
      if (actions != null) {
        assertEquals(1, actionsAlignmentScore(actions, false));
        assertEquals(0, actions[ActionsHelper.TEMPLATE_START_INDEX]);
        assertEquals("=========================NNNNN=========X", ActionsHelper.toString(actions));
      }

      actions = mCgGotoh.calculateEditDistance(s1, s1.length, t2, 0, 10, 5, true);
      if (actions != null) {
        assertEquals(5, actionsAlignmentScore(actions, false));
        assertEquals(0, actions[ActionsHelper.TEMPLATE_START_INDEX]);
        assertEquals("=========================NNNNN=====XXXXX", ActionsHelper.toString(actions));
      }
    } finally {
      Diagnostic.setLogStream();
    }
  }

  public void testalign3Mismatches() {
    getEditDistanceInstance(1, 1, 0);

    final byte[] r = DnaUtils.encodeString("gatcatatatCtatatatatatata     tatataaaac".replaceAll(" ", ""));
    final byte[] t = DnaUtils.encodeString("gatcatataCatatatatatCtata.....tatataaaac");

    int[] actions = mCgGotoh.calculateEditDistance(r, r.length, t, 0, 10, 3, true);
    if (actions != null) {
      assertEquals("=========XX=========X====NNNNN==========", ActionsHelper.toString(actions));
      assertEquals(3, actionsAlignmentScore(actions, false));
      assertEquals(0, actions[ActionsHelper.TEMPLATE_START_INDEX]);
    }
    getEditDistanceInstance(1, 1, 1);
    actions = mCgGotoh.calculateEditDistance(r, r.length, t, 0, 10, 3, true);
    if (actions != null) {
      assertEquals("=========XX=========X====NNNNN==========", ActionsHelper.toString(actions));
      assertEquals(3, actionsAlignmentScore(actions, false));
      assertEquals(0, actions[ActionsHelper.TEMPLATE_START_INDEX]);
    }
  }

  public void testAlignNs() {
    getEditDistanceInstance(1, 1, 0);

    final byte[] r = DnaUtils.encodeString("gatcatatatNtatatatatNtata     tatataaaac".replaceAll(" ", ""));
    final byte[] t = DnaUtils.encodeString("gatcatataNatatatatatNtata.....tatataaaac");

    int[] actions = mCgGotoh.calculateEditDistance(r, r.length, t, 0, 10, 3, true);
    if (actions != null) {
      assertEquals("=========TR=========R====NNNNN==========", ActionsHelper.toString(actions));
      assertEquals(0, actionsAlignmentScore(actions, false));
      assertEquals(0, actions[ActionsHelper.ALIGNMENT_SCORE_INDEX]);
      assertEquals(0, actions[ActionsHelper.TEMPLATE_START_INDEX]);
    }
    getEditDistanceInstance(1, 1, 1);
    actions = mCgGotoh.calculateEditDistance(r, r.length, t, 0, 10, 3, true);
    if (actions != null) {
      assertEquals("=========TR=========R====NNNNN==========", ActionsHelper.toString(actions));
      assertEquals(3, actionsAlignmentScore(actions, true));
      assertEquals(3, actions[ActionsHelper.ALIGNMENT_SCORE_INDEX]);
      assertEquals(0, actions[ActionsHelper.TEMPLATE_START_INDEX]);
    }
  }

  public void testAlignmentScores() {

    getEditDistanceInstance(1, 1, 0, "cg_test_errors");
//    Claimed alignment is incorrect:
//      ACCTCCTATGAAAAAACTTCCTACCACTCAC-CCTAG <-- template
//      ||||||||||      ||||||||||  ||| ||| |
//      ACCTCCTATG......CTTCCTACCAAACACTNCTTG <-- read
//      Observed mismatches: 4 Claimed mismatches: null
//      Observed alignment score: 5 Claimed alignment score: 3
//      Alignment mismatch 644988       137     chr1    554579  0       10=6N10=2X3=1I3=1X1=    *       0       0       ACCTCCTATGCTTCCTACCAAACACTNCTTG 6877-4.&04;:::;;7;:::89998!9988 AS:i:3
//    XS:Z:10=6N10=2X3=1I1X2I1=2B3=1X1=       XR:Z:AATGATT    XQ:Z:5578       IH:i:2  NH:i:2

    byte[] r = DnaUtils.encodeString("ACCTCCTATG     CTTCCTACCAAACACTGATCNCTTG".replaceAll(" ", ""));
    byte[] t = DnaUtils.encodeString("ACCTCCTATGAAAAAACTTCCTACCACTCACCCTAGGGGGGGGGGGG");
    int[] actions = mCgGotoh.calculateEditDistance(r, r.length, t, 0, 10, 3, false);

//    System.err.println(ActionsHelper.toString(actions) + " as=" + ActionsHelper.alignmentScore(actions));
    assertEquals(9, ActionsHelper.alignmentScore(actions));




//    GAGGAGGAGGGAATGTTTCCAAACTCAGATTTCTATGAGGCC <-- template
//    || || ||||| |   ||||||||||      ||||||||||
//    GATGNTGAGGGTAA..TTCCAAACTC......CTATGAGGCC <-- read
//    Observed mismatches: 4 Claimed mismatches: null
//    Observed alignment score: 4 Claimed alignment score: 3
//    Alignment mismatch 727966       153     chr22   45958796        37      2=1X2=1X5=1X1=1X2N10=6N10=      *       0       0 GATGNTGAGGGTAATTCCAAACTCCTATGAGGCC       899:!87989999::;;<;:98:9-0*2/74384      AS:i:3
//    XS:Z:2=1X2=2B1=1D1X5=1X1=1X2N10=6N10=      XR:Z:TTTA       XQ:Z:6  IH:i:1  NH:i:1
    r = DnaUtils.encodeString("GATGNGTGAGGGTAATTCCAAACTC     CTATGAGGCC".replaceAll(" ", ""));
    t = DnaUtils.encodeString("GAGGAGGAGGGAATGTTTCCAAACTCAGATTTCTATGAGGCC");

    actions = mCgGotoh.calculateEditDistance(r, r.length, t, 0, 10, 3, true);

//    System.err.println(ActionsHelper.toString(actions) + " as=" + ActionsHelper.alignmentScore(actions));
    assertEquals(6, ActionsHelper.alignmentScore(actions));


//    Claimed alignment is incorrect:
//      GAAGGTTGGGGGAGTGGGAGTTGGTGGCCTACCACTG-GGT <-- template
//      ||||||||||       |||||||||||||||||||| |||
//      GAAGGTTGNG.......GAGTTGGTNGCCTACCACTGTGGT <-- read
//      Observed alignment score: 4 Claimed alignment score: 3
//      Alignment mismatch 62705        131     paolo-bac       698     55      8=1X1=7N8=1X11=1I3=     paolo-bac       363     -335    GAAGGTTGNGGAGTTGGTNGCCTACCACTGTGGT      /21-/704!3.062/510!1,,,1042+/-5021      AS:i:3  NM:i:3  MQ:i:255        GS:Z:GG GC:Z:29S1G4S    GQ:Z:0  XA:i:3  IH:i:1  NH:i:1
    r = DnaUtils.encodeString("GAAGGTTGNG     GAGTTGGTNGCCTACCACTGGTGGT".replaceAll(" ", ""));
    t = DnaUtils.encodeString("GAAGGTTGGGGGAGTGGGAGTTGGTGGCCTACCACTGGGT");


    getEditDistanceInstance(1, 1, 1);
    actions = mCgGotoh.calculateEditDistance(r, r.length, t, 0, 10, 3, false);

//    System.err.println(ActionsHelper.toString(actions) + " as=" + ActionsHelper.alignmentScore(actions));
    assertEquals(4, ActionsHelper.alignmentScore(actions));

  }

  private int actionsAlignmentScore(int[] actions, boolean treatNsAsMismatches) {
    final ActionsHelper.CommandIterator it = ActionsHelper.iterator(actions);
    int prevAction = ActionsHelper.SAME;
    int ascore = 0;
    while (it.hasNext()) {
      final int currAction = it.next();
      if (currAction == ActionsHelper.MISMATCH) {
        ascore++;
      } else if (currAction == ActionsHelper.DELETION_FROM_REFERENCE) {
        if (prevAction != ActionsHelper.DELETION_FROM_REFERENCE) {
          ascore += 2;
        } else {
          ascore++;
        }
      } else if (currAction == ActionsHelper.INSERTION_INTO_REFERENCE) {
        if (prevAction != ActionsHelper.INSERTION_INTO_REFERENCE) {
          ascore += 2;
        } else {
          ascore++;
        }
      } else if (currAction == ActionsHelper.UNKNOWN_TEMPLATE || currAction == ActionsHelper.UNKNOWN_READ) {
        if (treatNsAsMismatches) {
          ascore++;
        }
      }
      prevAction = currAction;
    }
    return ascore;
  }
}

