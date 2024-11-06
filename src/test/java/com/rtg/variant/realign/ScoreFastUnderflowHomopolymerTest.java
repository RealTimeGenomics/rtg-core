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

import java.io.IOException;
import java.io.StringReader;
import java.util.Arrays;

import com.rtg.util.integrity.Exam;
import com.rtg.variant.util.arithmetic.LogApproximatePossibility;
import com.rtg.variant.util.arithmetic.SimplePossibility;

import junit.framework.TestCase;

/**
 */
public class ScoreFastUnderflowHomopolymerTest extends TestCase {

  private AllPaths score(Environment env) throws IOException {
    final AllPaths score = new ScoreFastUnderflowHomopolymer(new ScoreMatrixTest.MockRealignParams(), new HomoPolymerParams(SimplePossibility.SINGLETON, 3, 3, new StringReader("")),  new HomoPolymerParams(LogApproximatePossibility.SINGLETON, 3, 3, new StringReader("")));
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

  public void testUnderflow() throws IOException {
    final AllPaths score = score(envUnderflow());
    Exam.integrity(score);
    //System.err.println(score.totalScoreLn());
    //System.err.println(score.toString());
    assertEquals(-750.22187, score.totalScoreLn(), 0.001);
  }
}
