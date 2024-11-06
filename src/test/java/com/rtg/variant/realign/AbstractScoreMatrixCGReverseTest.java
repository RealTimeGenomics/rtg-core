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

import com.rtg.launcher.AbstractNanoTest;
import com.rtg.util.integrity.Exam;
import com.rtg.util.machine.MachineType;
import com.rtg.variant.util.arithmetic.LogApproximatePossibility;
import com.rtg.variant.util.arithmetic.LogPossibility;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;
import com.rtg.variant.util.arithmetic.SimplePossibility;


/**
 */
public abstract class AbstractScoreMatrixCGReverseTest extends AbstractNanoTest {

  //Mark says I have to put my name here so he doesnt get blamed for the next bit of code
  //JC

  /** Used for testing. */
  public static class LogTest extends AbstractScoreMatrixCGReverseTest {
    public LogTest() { }
    @Override
    protected PossibilityArithmetic arith() {
      return LogPossibility.SINGLETON;
    }
  }

  /** Used for testing. */
  public static class LogApproximateTest extends AbstractScoreMatrixCGReverseTest {
    public LogApproximateTest() { }
    @Override
    protected PossibilityArithmetic arith() {
      return LogApproximatePossibility.SINGLETON;
    }
    @Override
    public void testToString() { }
    @Override
    public void testToStringV2() { }
  }

  /** Used for testing. */
  public static class SimpleTest extends AbstractScoreMatrixCGReverseTest {
    public SimpleTest() { }
    @Override
    protected PossibilityArithmetic arith() {
      return SimplePossibility.SINGLETON;
    }
  }

  protected abstract PossibilityArithmetic arith();

  protected AllPaths score(final Environment env) {
    final AllPaths score = new ScoreMatrixCGReverse(arith(), new AbstractScoreMatrixCGTest.MockRealignParamsCG());
    score.setEnv(env);
    Exam.globalIntegrity(score);
    return score;
  }


  public void testRowOffset() {
    final AllPaths score0 = score(AbstractScoreMatrixCGTest.env());
    if (score0 instanceof ScoreMatrixCG) {
      final ScoreMatrixCG score = (ScoreMatrixCG) score0;
      final int half = -8; // half of the width of the band (at row 1).
      assertEquals(half, score.rowOffset(0));
      assertEquals(half + 5, score.rowOffset(5));
      assertEquals(half + 6 - 2, score.rowOffset(6));
      assertEquals(half + 25 - 2, score.rowOffset(25));
      assertEquals(half + 26 - 2 + 6, score.rowOffset(26));
    }
  }

  public void testStartScores() {
    final Environment env = AbstractScoreMatrixCGTest.env();
    final AllPaths score0 = score(env);
    if (score0 instanceof ScoreMatrixCG) {
      final ScoreMatrixCG score = (ScoreMatrixCG) score0;
      final int maxShift = env.maxShift();
      final int w = 2 * maxShift + 1;  // width
      assertEquals(w, score.mWidth);
      final int last = score.mLength;
      for (int col = 0; col < w; ++col) {
        assertEquals(0.0, score.arithmetic().poss2Ln(score.delete(last, col)));
        assertEquals(0.0, score.arithmetic().poss2Ln(score.match(last, col)));
        assertEquals(Double.NEGATIVE_INFINITY, score.arithmetic().poss2Ln(score.insert(last, col)));
      }
    }
  }

  public void testTotalScoreLn() {
    final AllPaths score = score(AbstractScoreMatrixCGTest.env());
    assertEquals(-10.012458, score.totalScoreLn(), 0.0002);
    assertEquals(-10.012458, Math.log(score.totalScore()), 0.0002);
  }

  public void testScoreLn() {
    final Environment env = AbstractScoreMatrixCGTest.env();
    final AllPaths score0 = score(env);
    if (score0 instanceof ScoreMatrixCG) {
      final ScoreMatrixCG score = (ScoreMatrixCG) score0;
      final int maxShift = env.maxShift();
      assertEquals(-maxShift, score.rowOffset(1));
      final int start = env.absoluteTemplatePosition(0);
      // down column 8 (first T on the template), cuts through the overlap region
      final int tPos = 8 - start;
      AbstractScoreMatrixCGTest.checkScores(score, 0, tPos, -12.6462, -12.4498, Double.NEGATIVE_INFINITY);
      AbstractScoreMatrixCGTest.checkScores(score, 1, tPos, -16.8551, -16.9501, -13.0428);
      AbstractScoreMatrixCGTest.checkScores(score, 2, tPos, -14.7730, -14.4005, -10.3354);
      AbstractScoreMatrixCGTest.checkScores(score, 3, tPos, -16.5896, -16.1711, -11.9198);
      AbstractScoreMatrixCGTest.checkScores(score, 4, tPos, -19.0452, -18.3808, -13.5117);
      AbstractScoreMatrixCGTest.checkScores(score, 5, tPos, -21.3363, -20.4260, -15.0958);
      AbstractScoreMatrixCGTest.checkScores(score, 6, tPos, -19.7458, -18.8543, -13.7007);
      AbstractScoreMatrixCGTest.checkScores(score, 7, tPos, -26.3596, -21.1361, -15.2996);
      AbstractScoreMatrixCGTest.checkScores(score, 8, tPos, -27.9488, -22.7212, -16.8846);
      AbstractScoreMatrixCGTest.checkScores(score, 9, tPos, -29.5036, -24.3107, -18.4742);
      AbstractScoreMatrixCGTest.checkScores(score, 10, tPos, -29.8945, -25.8869, -20.0569);
      AbstractScoreMatrixCGTest.checkScores(score, 11, tPos, -27.6707, -26.7838, -21.6349);
    }
  }

  public void testToString() throws IOException {
    final AllPaths score = score(AbstractScoreMatrixCGTest.env());
    mNano.check("scorematrixcg-reverse-tostring.txt", score.toString());
  }

  protected AllPaths scoreV2(final Environment env) {
    final AllPaths score = new ScoreMatrixCG(arith(), new AbstractScoreMatrixCGTest.MockRealignParamsCG(MachineType.COMPLETE_GENOMICS_2));
    score.setEnv(env);
    Exam.globalIntegrity(score);
    return score;
  }

  public void testToStringV2() throws IOException {
    final AllPaths score = scoreV2(AbstractScoreMatrixCGTest.envV2());
    mNano.check("scorematrixcg-reverse-v2-tostring.txt", score.toString());
  }

  public void testReadStartsBeforeLn() {
    final AllPaths score = score(AbstractScoreMatrixCGTest.env());
    if (score instanceof ScoreMatrixCG) {
      final ScoreMatrixCGReverse sm = (ScoreMatrixCGReverse) score;
      assertEquals(Double.NEGATIVE_INFINITY, sm.readStartsBeforeLn(-1));
      assertEquals(-23.9712, sm.readStartsBeforeLn(0), 0.0001);
      assertEquals(-22.1798, sm.readStartsBeforeLn(1), 0.0001);
      assertEquals(-20.5382, sm.readStartsBeforeLn(2), 0.0001);
      assertEquals(-18.9231, sm.readStartsBeforeLn(3), 0.0001);
      assertEquals(-17.3115, sm.readStartsBeforeLn(4), 0.0001);
      assertEquals(-15.4648, sm.readStartsBeforeLn(5), 0.0001);
      assertEquals(-14.0367, sm.readStartsBeforeLn(6), 0.0001);
      assertEquals(-4.5470, sm.readStartsBeforeLn(7), 0.0001);
      assertEquals(-4.5435, sm.readStartsBeforeLn(8), 0.0001);
      assertEquals(-2.3226, sm.readStartsBeforeLn(9), 0.0001);
      assertEquals(-0.0140, sm.readStartsBeforeLn(10), 0.0001);
      assertEquals(-0.0013, sm.readStartsBeforeLn(11), 0.0001);
      assertEquals(-0.0001, sm.readStartsBeforeLn(12), 0.0001);
      assertEquals(0.0000, sm.readStartsBeforeLn(13), 0.0001);
      assertEquals(0.0000, sm.readStartsBeforeLn(14), 0.0001);
      assertEquals(0.0, sm.readStartsBeforeLn(15));
    }
  }
}
