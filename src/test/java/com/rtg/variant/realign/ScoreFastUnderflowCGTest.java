/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

package com.rtg.variant.realign;

import java.util.Arrays;

import com.rtg.util.integrity.Exam;
import com.rtg.util.machine.MachineType;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;

/**
 */
public class ScoreFastUnderflowCGTest extends AbstractScoreMatrixCGTest {

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
