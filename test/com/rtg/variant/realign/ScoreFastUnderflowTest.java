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

import java.util.Arrays;

import com.rtg.util.integrity.Exam;

import junit.framework.Test;
import junit.framework.TestSuite;

/**
 */
public class ScoreFastUnderflowTest extends ScoreMatrixTest {

  public static Test suite() {
    final TestSuite suite = new TestSuite();
    suite.addTestSuite(ScoreFastUnderflowTest.class);
    return suite;
  }

  @Override
  protected AllPaths score(Environment env) {
    final AllPaths score = new ScoreFastUnderflow(new MockRealignParams());
    score.setEnv(env);
    Exam.globalIntegrity(score);
    return score;
  }

  protected Environment envUnderflow() {
    //an environment that should cause an underflow
    final int length = 150;
    final byte[] read = new byte[length];
    Arrays.fill(read, (byte) 1);
    final byte[] template = new byte[2 * length];
    Arrays.fill(template, (byte) 2);
    final int start = 2;
    final int maxShift = 2;
    final double[] quality = new double[length];
    Arrays.fill(quality, 0.01);
    final Environment env = new EnvironmentImplementation(
        maxShift,
        template,
        start,
        read,
        quality
        );
    Exam.integrity(env);
    return env;
  }

  public void testUnderflow() {
    final AllPaths score = score(envUnderflow());
    Exam.integrity(score);
    //System.err.println(score.totalScoreLn());
    //System.err.println(score.toString());
    final double expected = -750.22187;
    assertEquals(expected, score.totalScoreLn(), 0.001);
    //assertEquals(expected, Math.log(score.totalScore()), 0.001); //underflows
    assertEquals(expected, score.arithmetic().poss2Ln(score.total()), 0.001);
  }
}
