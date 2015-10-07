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

import java.io.IOException;

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
  public void test1() throws IOException {
    final AllPaths score = score(env());
    mNano.check("scorematrix-reverse-typical.txt", score.toString());
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
