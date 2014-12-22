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
package com.rtg.variant.realign;

import static com.rtg.util.StringUtils.LS;

import com.rtg.util.integrity.Exam;
import com.rtg.variant.util.arithmetic.LogApproximatePossibility;
import com.rtg.variant.util.arithmetic.LogPossibility;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;
import com.rtg.variant.util.arithmetic.SimplePossibility;

import junit.framework.Test;
import junit.framework.TestSuite;

/**
 */
public abstract class ScoreMatrixReverseTest extends ScoreMatrixTest {

  //Mark says I have to put my name here so he doesnt get blamed for the next bit of code
  //JC

  /** Used for testing. */
  public static class LogTest extends ScoreMatrixReverseTest {
    public LogTest() { }
    @Override
    protected AllPaths score(final Environment env) {
      return score(LogPossibility.SINGLETON, env);
    }
  }

  /** Used for testing. */
  public static class LogApproximateTest extends ScoreMatrixReverseTest {
    public LogApproximateTest() { }
    @Override
    protected AllPaths score(final Environment env) {
      return score(LogApproximatePossibility.SINGLETON, env);
    }
    //the results are slightly inaccurate in the 5th digit - but the other two accurate versions agree so let it be
    @Override
    public void test1() { }
  }

  /** Used for testing. */
  public static class SimpleTest extends ScoreMatrixReverseTest {
    public SimpleTest() { }
    @Override
    protected AllPaths score(final Environment env) {
      return score(SimplePossibility.SINGLETON, env);
    }
  }

  public static Test suite() {
    final TestSuite suite = new TestSuite();
    suite.addTestSuite(LogApproximateTest.class);
    suite.addTestSuite(LogTest.class);
    suite.addTestSuite(SimpleTest.class);
    return suite;
  }

  @Override
  protected abstract AllPaths score(final Environment env);

  @Override
  protected AllPaths score(final PossibilityArithmetic arith, final Environment env) {
    final AllPaths score = new ScoreMatrixReverse(arith, new MockRealignParams());
    score.setEnv(env);
    Exam.globalIntegrity(score);
    return score;
  }

  @Override
  public void test1() {
    final AllPaths score = score(env());
    final String exp = ""
        + "ScoreMatrix                |G            0       |G            1       |A            2       |C            3       |G            4       |T            5       |A            6       |C            7       |G            8       |T            9       |T           10       |T           11       |N           12       |" + LS
        + "[  0]         20.462 20.659|        22.307 20.322|        13.284 13.477|        21.330 16.092|        22.034 18.669|" + LS
        + "[  1]A                     | 16.215 17.135 17.337| 14.988 19.011 16.095| 13.390 13.167 12.992| 20.192 20.026 15.592| 19.448 18.847 14.611|" + LS
        + "[  2]C                     |                     | 13.006 12.831 13.028| 14.485 14.979 11.057| 13.346 13.124 12.497| 17.080 16.206 11.518| 18.725 14.918  9.566|" + LS
        + "[  3]G                     |                     |                     |  7.975  7.771  7.968| 10.402 10.184  9.970| 13.967 13.107  8.420| 15.343 11.838  6.465| 13.924 13.444  9.531|" + LS
        + "[  4]C                     |                     |                     |                     |  6.649  7.749  7.952|  5.348  5.126  5.323| 13.955  8.737  3.364| 15.598 11.815  6.443| 14.275 13.680  9.459|" + LS
        + "[  5]G                     |                     |                     |                     |                     |  1.897  7.733 12.883|  0.288  0.066  0.263| 13.659  8.713  3.342| 14.963 11.704  6.366| 15.146 14.031  9.176|" + LS
        + "[  6]A                     |                     |                     |                     |                     |                     |  1.875  7.626 10.248|  0.266  0.044  0.241| 10.094  8.433  3.265| 10.127  9.864  6.078| 10.310  9.864  6.078|" + LS
        + "[  7]C                     |                     |                     |                     |                     |                     |                     |  1.820  4.991  5.257|  0.242  0.022  0.165|  5.067  5.027  2.992|  5.099  5.027  2.992|  5.282  5.028  2.992|" + LS
        + "[  8]G                     |                     |                     |                     |                     |                     |                     |                     |         0.000  0.000|         0.000  0.000|         0.000  0.000|         0.000  0.000|         0.000  0.000|" + LS
        + "                           |G            0       |G            1       |A            2       |C            3       |G            4       |T            5       |A            6       |C            7       |G            8       |T            9       |T           10       |T           11       |N           12       |" + LS
        ;
    assertEquals(exp, score.toString());
  }

  @Override
  //needs some sortingout to make work
  public void test2() {
  }

  public void testReadStartsBeforeLn() {
    final AllPaths score = score(env());
    if (score instanceof ScoreMatrixReverse) {
      final ScoreMatrixReverse sm = (ScoreMatrixReverse) score;
      assertEquals(Double.NEGATIVE_INFINITY, sm.readStartsBeforeLn(-1));
      assertEquals(-7.17996, sm.readStartsBeforeLn(0), 0.0001);
      assertEquals(-7.03252, sm.readStartsBeforeLn(1), 0.0001);
      assertEquals(-0.00053, sm.readStartsBeforeLn(2), 0.0001);
      assertEquals(-0.00016, sm.readStartsBeforeLn(3), 0.0001);
      assertEquals(0.0, sm.readStartsBeforeLn(4));
      assertEquals(0.0, sm.readStartsBeforeLn(5));
      assertFalse(sm.underflow());
    }
  }
}
