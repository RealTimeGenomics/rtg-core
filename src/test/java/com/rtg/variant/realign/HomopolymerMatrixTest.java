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

import static com.rtg.util.StringUtils.LS;

import java.io.IOException;
import java.io.Reader;
import java.io.StringReader;
import java.util.Arrays;

import com.rtg.launcher.AbstractNanoTest;
import com.rtg.mode.DnaUtils;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;
import com.rtg.util.machine.MachineType;
import com.rtg.variant.util.arithmetic.LogPossibility;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;

/**
 */
public class HomopolymerMatrixTest extends AbstractNanoTest {

  private static final String IONT = ""
      + "A+T 3 0 10 10 10" + LS
      + "A+T 4 0 0 10 10 10" + LS
      + "C+G 3 0 42 42 0" + LS
      + "C+G 4 0 0 0 42 42" + LS
      ;
  protected AllPaths score(final EnvironmentHomopolymer env) throws IOException {
    return score(LogPossibility.SINGLETON, env, IONT);
  }

  protected AllPaths score(final PossibilityArithmetic arith, final EnvironmentHomopolymer env, final String calibration) throws IOException {
    final Reader in = new StringReader(calibration.replaceAll(" ", "\t"));
    final HomoPolymerParams params = new HomoPolymerParams(arith, 2, 3, in);
    //System.err.println(params);
    //System.err.println(params.transition(2, 3, 2));
    final AllPaths score = new HomopolymerMatrix(arith, new MockRealignParams(), params);
    score.setEnv(env);
    Exam.globalIntegrity(score);
    return score;
  }

  /**
   * Realign parameters that match the scorematrixtestCG.xls spreadsheet.
   */
  public static class MockRealignParams extends IntegralAbstract implements RealignParams {
    @Override
    public double matchLn() {
      return Math.log(0.98939);
    }
    @Override
    public double misMatchLn() {
      return Math.log(0.00920);
    }
    @Override
    public double insertOpenLn() {
      return Math.log(0.00058);
    }
    @Override
    public double insertExtendLn() {
      return Math.log(0.2);
    }
    @Override
    public double deleteOpenLn() {
      return Math.log(0.00083);
    }
    @Override
    public double deleteExtendLn() {
      return Math.log(0.18);
    }
    @Override
    public MachineType machineType() {
      return null;
    }
    @Override
    public int gapEnd(final int gap) {
      throw new UnsupportedOperationException();
    }
    @Override
    public double gapFreqLn(final int gap, final int width) {
      throw new UnsupportedOperationException();
    }
    @Override
    public int gapStart(final int gap) {
      throw new UnsupportedOperationException();
    }
    @Override
    public void toString(final StringBuilder sb) {
      // do nothing
    }
    @Override
    public boolean integrity() {
      return true;
    }
    @Override
    public double[][] gapDistributionPoss(PossibilityArithmetic arith) {
      throw new UnsupportedOperationException();
    }
  }

  protected EnvironmentHomopolymer env(final String read, final String template) {
    //see scorematrixtest.xls for details
    final int start = 2;
    final int maxShift = 2;
    final byte[] readBytes = DnaUtils.encodeString(read);
    final double[] quality = new double[readBytes.length];
    Arrays.fill(quality, 0.01);
    quality[0] = 0.1;
    quality[1] = Math.pow(10.0, -15 / 10.0);
    final byte[] templateBytes = DnaUtils.encodeString(template);
    final Environment env = new EnvironmentImplementation(
        maxShift,
        templateBytes,
        start,
        readBytes,
        quality
        );
    Exam.integrity(env);
    return new EnvironmentHomopolymer(env);
  }


  public void test1() throws IOException {
    final EnvironmentHomopolymer env = env(ScoreMatrixTest.READ, ScoreMatrixTest.TEMPLATE);
    //System.err.println(env);
    final AllPaths score = score(env);
    mNano.check("homopolymermatrix-test.txt", score.toString());
  }

  //homopolymer reads that invoke short circuit calculations in HomopolymerMatrix
  public void test2a() throws IOException {
    final EnvironmentHomopolymer env = env(ScoreMatrixTest.READ_H, ScoreMatrixTest.TEMPLATE_H);
    //System.err.println(env);
    final AllPaths score = score(LogPossibility.SINGLETON, env, IONT);
    mNano.check("homopolymermatrix-short-circuit.txt", score.toString());
  }

  //homopolymer reads that invoke short circuit calculations in HomopolymerMatrix but with calibration table empty so should be same as ScoreMatrix results
  public void test2b() throws IOException {
    final EnvironmentHomopolymer env = env(ScoreMatrixTest.READ_H, ScoreMatrixTest.TEMPLATE_H);
    //System.err.println(env);
    final AllPaths score = score(LogPossibility.SINGLETON, env, "");
    mNano.check("scorematrix-homopolymer.txt", score.toString()); // This expected result is shared with ScoreMatrix
  }

  public void testTotalScoreLn() throws IOException {
    final EnvironmentHomopolymer env = env(ScoreMatrixTest.READ_H, ScoreMatrixTest.TEMPLATE_H);
    final AllPaths score = score(LogPossibility.SINGLETON, env, IONT);
    //regression tests
    assertEquals(-4.44028, score.totalScoreLn(), 0.0001);
    assertEquals(-4.44028, Math.log(score.totalScore()), 0.0001);
    if (score instanceof AbstractAllPaths) {
      assertEquals(-env.maxShift(), ((AbstractAllPaths) score).rowOffset(1));
    }
  }

  public void testReadEndsAfterLn() throws IOException {
    final AllPaths score = score(LogPossibility.SINGLETON, env(ScoreMatrixTest.READ_H, ScoreMatrixTest.TEMPLATE_H), IONT);
    if (score instanceof ScoreMatrix) {
      final ScoreMatrix sm = (ScoreMatrix) score;
      //regression tests
      assertEquals(0.0, sm.readEndsAfterLn(-1));
      assertEquals(0.0, sm.readEndsAfterLn(7));
      assertEquals(-7.3387e-7, sm.readEndsAfterLn(8), 0.0001e-7);
      assertEquals(-1.5639e-5, sm.readEndsAfterLn(9), 0.0001e-5);
      assertEquals(-4.9783e-4, sm.readEndsAfterLn(10), 0.0001e-4);
      assertEquals(-5.7572, sm.readEndsAfterLn(11), 0.0001);
      assertEquals(Double.NEGATIVE_INFINITY, sm.readEndsAfterLn(12));
      assertFalse(sm.underflow());
    }
  }

}
