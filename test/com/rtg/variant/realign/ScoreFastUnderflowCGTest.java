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
import com.rtg.util.machine.MachineType;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;

/**
 */
public class ScoreFastUnderflowCGTest extends ScoreMatrixCGTest {

  @Override
  protected PossibilityArithmetic arith() {
    return null;
  }

  @Override
  protected AllPaths score(Environment env) {
    final AllPaths score = new ScoreFastUnderflowCG(new MockRealignParamsCG());
    score.setEnv(env);
    Exam.globalIntegrity(score);
    return score;
  }

  @Override
  protected AllPaths scoreV2(final Environment env) {
    final AllPaths score = new ScoreFastUnderflowCG(new MockRealignParamsCG(MachineType.COMPLETE_GENOMICS_2));
    score.setEnv(env);
    Exam.globalIntegrity(score);
    return score;
  }

  protected Environment envUnderflow() {
    //an environment that should cause an underflow
    final int length = 35;
    final byte[] read = new byte[length];
    Arrays.fill(read, (byte) 1);
    final byte[] template = new byte[2 * length];
    Arrays.fill(template, (byte) 2);
    final int start = 2;
    final int maxShift = 2;
    final double[] quality = new double[length];
    Arrays.fill(quality, 0.000001);
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
    //it does take a bit of work but it is possible to get cg to underflow
    final MockRealignParamsCG params = new MockRealignParamsCG() {
      @Override
      public double deleteOpenLn() {
        return -1000.0;
      }
      @Override
      public double insertOpenLn() {
        return -1000.0;
      }
      @Override
      public double matchLn() {
        return -1000.0;
      }
      @Override
      public double misMatchLn() {
        return -1000.0;
      }
    };
    final AllPaths score1 = new ScoreFastUnderflowCG(params);
    score1.setEnv(envUnderflow());
    //System.err.println(score1.toString());
    Exam.globalIntegrity(score1);
    Exam.integrity(score1);
    //System.err.println(score.totalScoreLn());
    //System.err.println(score.toString());
    assertEquals(-22101.3, score1.totalScoreLn(), 0.1);
  }
}
