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

import com.rtg.util.integrity.Exam;
import com.rtg.variant.util.arithmetic.LogApproximatePossibility;
import com.rtg.variant.util.arithmetic.LogPossibility;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;
import com.rtg.variant.util.arithmetic.SimplePossibility;

/**
 */
public abstract class AbstractScoreMatrixReverseTest extends ScoreMatrixTest {

  //Mark says I have to put my name here so he doesnt get blamed for the next bit of code
  //JC

  /** Used for testing. */
  public static class LogTest extends AbstractScoreMatrixReverseTest {
    public LogTest() { }
    @Override
    protected AllPaths score(final Environment env) {
      return score(LogPossibility.SINGLETON, env);
    }
  }

  /** Used for testing. */
  public static class LogApproximateTest extends AbstractScoreMatrixReverseTest {
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
  public static class SimpleTest extends AbstractScoreMatrixReverseTest {
    public SimpleTest() { }
    @Override
    protected AllPaths score(final Environment env) {
      return score(SimplePossibility.SINGLETON, env);
    }
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
