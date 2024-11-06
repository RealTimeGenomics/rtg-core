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
import java.util.Arrays;

import com.rtg.launcher.AbstractNanoTest;
import com.rtg.mode.DnaUtils;
import com.rtg.util.integrity.Exam;
import com.rtg.util.machine.MachineType;
import com.rtg.variant.util.arithmetic.LogApproximatePossibility;
import com.rtg.variant.util.arithmetic.LogPossibility;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;
import com.rtg.variant.util.arithmetic.SimplePossibility;

/**
 */
public abstract class AbstractScoreMatrixCGTest extends AbstractNanoTest {

  //Mark says I have to put my name here so he doesnt get blamed for the next bit of code
  //JC
  /** Used for testing. */
  public static class LogTest extends AbstractScoreMatrixCGTest {
    public LogTest() { }
    @Override
    protected PossibilityArithmetic arith() {
      return LogPossibility.SINGLETON;
    }
  }

  /** Used for testing. */
  public static class LogApproximateTest extends AbstractScoreMatrixCGTest {
    public LogApproximateTest() { }
    @Override
    protected PossibilityArithmetic arith() {
      return LogApproximatePossibility.SINGLETON;
    }
    //the results are slightly inaccurate in the 5th digit - but the other two accurate versions agree so let it be
    @Override
    public void testToString() { }
    @Override
    public void testToStringV2() { }
  }

  /** Used for testing. */
  public static class SimpleTest extends AbstractScoreMatrixCGTest {
    public SimpleTest() { }
    @Override
    protected PossibilityArithmetic arith() {
      return SimplePossibility.SINGLETON;
    }
  }

  /**
   * CG Realign parameters that match the scorematrixtestCG.xls spreadsheet.
   */
  public static class MockRealignParamsCG extends ScoreMatrixTest.MockRealignParams {

    static final int[] GAP_START = {-3, 0, 5, -7};

    static final double[][] GAP_FREQ =
      {
      {Math.log(0.08), Math.log(0.84), Math.log(0.08)},    // for gap of -3, -2, -1
      {Math.log(0.27), Math.log(0.635), Math.log(0.095)},  // for gap of  0,  1,  2
      {Math.log(0.90), Math.log(0.07), Math.log(0.03)},    // for gap of  5,  6,  7
      {Math.log(0.025), Math.log(0.047), Math.log(0.104), Math.log(0.306), Math.log(0.445), Math.log(0.07), Math.log(0.003)},    // for gap of -7, -6, -5, -4, -3, -2, -1
      };

    final MachineType mMachineType;

    public MockRealignParamsCG() {
      this(MachineType.COMPLETE_GENOMICS);
    }

    MockRealignParamsCG(MachineType mt) {
      mMachineType = mt;
    }

    @Override
    public MachineType machineType() {
      return mMachineType;
    }

    @Override
    public int gapStart(final int gap) {
      return GAP_START[gap];
    }

    @Override
    public int gapEnd(final int gap) {
      return GAP_START[gap] + GAP_FREQ[gap].length - 1;
    }

    @Override
    public double gapFreqLn(final int gap, final int width) {
      return GAP_FREQ[gap][width - gapStart(gap)];
    }

    @Override
    public double[][] gapDistributionPoss(PossibilityArithmetic arith) {
      return RealignParamsImplementation.gapDistributionPoss(GAP_FREQ, arith);
    }
  }

  /**
   * This is the template and read from the <code>scorematrixtestCG.xls</code> spreadsheet.
   * The read has an deletion of 3xA, with fragments having overlap of 2, a small gap of 0, and a large gap of 7.
   */
  public static final String TEMPLATE = "GGGGGGGATA AAAAA   GGCGACAT GCCAATGTGT CGCCTTT TTCAACTTTC CGATTAA".replaceAll(" ", "");
  static final String READ = "                  ATAAA     AAGGCGACAT GCCAATGTGT         TTCAACTTTC".replaceAll(" ", "");

  /*
                                              \del/
                                      GGGGGGGATAAAAAAGGCGACATGCCAATGTGTCGCCTTTTTCAACTTTCCGATTAA
Possible alignments:
Good overlap, and "actual" deletion of 3xA
                                             AT---AAA
                                                   AAGGCGACAT
                                                             GCCAATGTGT
                                                                              TTCAACTTTC
No mismatch, but unrealistic overlap
                                             ATAAA
                                                   AAGGCGACAT
                                                             GCCAATGTGT
                                                                              TTCAACTTTC
Likely overlap, 1 Substitution
                                                ATAAA
                                                   AAGGCGACAT
                                                             GCCAATGTGT
                                                                              TTCAACTTTC

29bp equivalent case
                                              \del/
                                      GGGGGGGATAAAAAAAAAAAGGCGACATGCCAATGTGTCGCCT
Good overlap, and "actual" deletion of 3xA
                                             AT---AAAAAAAA
                                                       AAAGGCGAC
                                                                ATGCCAATGT
No mismatch, but unrealistic overlap
                                             ATAAAAAAAA
                                                       AAAGGCGAC
                                                                ATGCCAATGT
Likely overlap, 1 Substitution
                                                ATAAAAAAAA
                                                       AAAGGCGAC
                                                                ATGCCAATGT
  */

  /**
   * This is the 29 bp equivalent case.
   * The read has a deletion of 3xA, with fragments having overlap of 3, and a small gap of 0.
   */
  public static final String TEMPLATE_2 = "GGGGGGGAT AAAAAAAAAAAGGCGAC ATGCCAATGT GTCGCCT".replaceAll(" ", "");
  static final String READ_2 = "                  ATAAAAAAAA AAAGGCGAC ATGCCAATGT".replaceAll(" ", "");


  protected static Environment env() {
    return env(DnaUtils.encodeString(READ), DnaUtils.encodeString(TEMPLATE));
  }

  protected static Environment env(byte[] read, byte[] template) {
    return env(read, template, 7);
  }

  protected static Environment env(byte[] read, byte[] template, int maxShift) {
    final double[] quality = new double[read.length];
    Arrays.fill(quality, 0.01);
    quality[0] = 0.1;
    quality[1] = Math.pow(10.0, -15 / 10.0);
    final Environment env = new EnvironmentImplementation(
        maxShift,
        template,
        7, //start
        read,
        quality
        );
    Exam.integrity(env);
    return env;
  }

  protected abstract PossibilityArithmetic arith();

  protected AllPaths score(final Environment env) {
    final AllPaths score = new ScoreMatrixCG(arith(), new MockRealignParamsCG());
    score.setEnv(env);
    Exam.globalIntegrity(score);
    return score;
  }

  protected static Environment envV2() {
    return env(DnaUtils.encodeString(READ_2), DnaUtils.encodeString(TEMPLATE_2), 9);
  }
  protected AllPaths scoreV2(final Environment env) {
    final AllPaths score = new ScoreMatrixCG(arith(), new MockRealignParamsCG(MachineType.COMPLETE_GENOMICS_2));
    score.setEnv(env);
    Exam.globalIntegrity(score);
    return score;
  }

  public void testRowOffset() {
    final AllPaths score0 = score(env());
    if (score0 instanceof ScoreMatrixCG) {
      final ScoreMatrixCG score = (ScoreMatrixCG) score0;
      final int half = -8; // half of the width of the band (at row 1).
      assertEquals(half, score.rowOffset(0));
      assertEquals(half + 1, score.rowOffset(1));
      assertEquals(half + 5, score.rowOffset(5));
      assertEquals(half + 6 - 2, score.rowOffset(6));
      assertEquals(half + 15 - 2, score.rowOffset(15));
      assertEquals(half + 16 - 2, score.rowOffset(16));
      assertEquals(half + 25 - 2, score.rowOffset(25));
      assertEquals(half + 26 - 2 + 6, score.rowOffset(26));
      assertEquals(half + 35 - 2 + 6 , score.rowOffset(35));
    }
  }

  public void testStartScores() {
    final Environment env = env();
    final AllPaths score0 = score(env);
    if (score0 instanceof ScoreMatrixCG) {
      final ScoreMatrixCG score = (ScoreMatrixCG) score0;
      final int maxShift = env.maxShift();
      final int w = 2 * maxShift + 1;  // width
      final double matchPenalty = Math.log(1.0 - Math.exp(score.mParams.deleteOpenLn())); // / w);
      final double deletePenalty = Math.log(Math.exp(score.mParams.deleteOpenLn())); // / w);
      assertEquals(w, score.mWidth);
      assertEquals(-maxShift - 1, score.rowOffset(0));
      for (int tPos = -maxShift - 1; tPos < maxShift; ++tPos) {
        checkScores(score, 0, tPos, deletePenalty, matchPenalty, Double.NEGATIVE_INFINITY);
      }
    }
  }

  public void testTotalScoreLn() {
    final AllPaths score = score(env());
    assertEquals(-10.012470, score.totalScoreLn(), 0.001);
    assertEquals(-10.012470, Math.log(score.totalScore()), 0.001);
  }

  public void testScoreLn() {
    final Environment env = env();
    final AllPaths score0 = score(env);
    if (score0 instanceof ScoreMatrixCG) {
      final ScoreMatrixCG score = (ScoreMatrixCG) score0;
      final int start = env.absoluteTemplatePosition(0);
      // down column 32 (T on the template), cuts through the big gap
      final int tPos = 32 - start;
      checkScores(score, 21, tPos, Double.NEGATIVE_INFINITY, -18.7048, -18.4873);
      checkScores(score, 22, tPos, -27.1852, -22.1531, -16.8977);
      checkScores(score, 23, tPos, -29.7517, -14.4359, -15.3150);
      checkScores(score, 24, tPos, -22.9162, -18.9807, -13.7276);
      checkScores(score, 25, tPos, -25.8054, -6.2972, -25.7807);
      checkScores(score, 26, tPos, -35.6482, -30.2646, -44.9473);
      checkScores(score, 27, tPos, -38.0540, -37.9438, -48.9078);
      checkScores(score, 28, tPos, -41.1499, -46.0732, -56.0204);
      checkScores(score, 29, tPos, -44.2510, -50.1387, Double.NEGATIVE_INFINITY);
    }
  }

  protected static void checkScores(ScoreMatrixCG score, final int row, final int relTemplatePos, final double delete, final double match, final double insert) {
    final int col = relTemplatePos - score.rowOffset(row);
    assertTrue("template position too small in test", 0 <= col);
    assertTrue("template position too large in test", col < score.mWidth);
    assertEquals(delete, score.arithmetic().poss2Ln(score.delete(row, col)), 0.001);
    assertEquals(match, score.arithmetic().poss2Ln(score.match(row, col)), 0.001);
    assertEquals(insert, score.arithmetic().poss2Ln(score.insert(row, col)), 0.001);
  }

  public void testToString() throws IOException {
    final AllPaths score = score(env());
    mNano.check("scorematrixcg-tostring.txt", score.toString());
  }

  public void testToStringV2() throws IOException {
    final AllPaths score = scoreV2(envV2());
    mNano.check("scorematrixcg-v2-tostring.txt", score.toString());
  }

  public void testReadEndsAfterLn() {
    final AllPaths score = score(env());
    if (score instanceof ScoreMatrixCG) {
      final ScoreMatrixCG sm = (ScoreMatrixCG) score;
      assertEquals(0.0, sm.readEndsAfterLn(37));
      assertEquals(0.0, sm.readEndsAfterLn(48), 0.0001);
      assertEquals(-0.0002, sm.readEndsAfterLn(49), 0.0001);
      assertEquals(-7.6540, sm.readEndsAfterLn(50), 0.0001);
      assertEquals(-14.1332, sm.readEndsAfterLn(51), 0.0001);
      assertEquals(-15.9250, sm.readEndsAfterLn(52), 0.0001);
      assertEquals(Double.NEGATIVE_INFINITY, sm.readEndsAfterLn(53));
    }
  }
}
