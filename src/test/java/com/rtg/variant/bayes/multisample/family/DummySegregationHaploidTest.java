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

package com.rtg.variant.bayes.multisample.family;

import com.rtg.util.MathUtils;
import com.rtg.util.integrity.Exam;
import com.rtg.variant.bayes.Code;
import com.rtg.variant.bayes.CodeDiploid;

import junit.framework.TestCase;

/**
 */
public class DummySegregationHaploidTest extends TestCase {

  public void testTrivial() {
    final Code code = new CodeDiploid(2);
    final ISegregationScore score = AbstractSegregationHaploid.getHaploidInstance(code, code.code(0), code.code(1, 1));
    assertTrue(score instanceof SegregationTrivial);
    Exam.integrity(score);
    score.increment(42, false);
    Exam.integrity(score);
    assertEquals(0.0, score.lnProbability());
  }

  public void testXY() {
    checkXY(0, 1);
    checkXY(1, 0);
  }

  private void checkXY(int x, int y) {
    final Code code = new CodeDiploid(3);
    final int q = 2;
    final ISegregationScore score = AbstractSegregationHaploid.getHaploidInstance(code, code.code(y), code.code(x, y));
    Exam.integrity(score);
    score.increment(code.code(x, y), true);
    Exam.integrity(score);
    score.increment(code.code(y, x), true);
    Exam.integrity(score);
    score.increment(code.code(y, y), true);
    Exam.integrity(score);
    score.increment(code.code(y, y), true);
    Exam.integrity(score);
    score.increment(code.code(y, y), true);
    Exam.integrity(score);
    score.increment(code.code(x), false);
    Exam.integrity(score);
    score.increment(code.code(y), false);
    Exam.integrity(score);

    //Should have no effect
    for (int i = 0; i < 10; ++i) {
      score.increment(code.code(x, x), true);
      Exam.integrity(score);
      score.increment(code.code(y, q), true);
      Exam.integrity(score);
      score.increment(code.code(q), false);
      Exam.integrity(score);
    }

    final double expected = MathUtils.logFactorial(7) - MathUtils.LOG_2 * 7 - MathUtils.logFactorial(3) - MathUtils.logFactorial(4);
    //System.err.println(expected);
    assertEquals(expected, score.lnProbability(), 0.00000001);
  }

  public void testXYZ() {
    checkXYX(0, 1, 2);
    checkXYX(0, 2, 1);
    checkXYX(1, 0, 2);
    checkXYX(1, 2, 0);
    checkXYX(2, 0, 1);
    checkXYX(2, 1, 0);
  }

  private void checkXYX(int x, int y, int z) {
    final Code code = new CodeDiploid(4);
    final int q = 3;
    final ISegregationScore score = AbstractSegregationHaploid.getHaploidInstance(code, code.code(z), code.code(x, y));
    Exam.integrity(score);
    score.increment(code.code(x, z), true);
    Exam.integrity(score);
    score.increment(code.code(y, z), true);
    Exam.integrity(score);
    score.increment(code.code(z, y), true);
    Exam.integrity(score);
    score.increment(code.code(x), false);
    Exam.integrity(score);
    score.increment(code.code(y), false);
    Exam.integrity(score);

    //Should have no effect
    for (int i = 0; i < 10; ++i) {
      score.increment(code.code(x, y), true);
      Exam.integrity(score);
      score.increment(code.code(y, x), true);
      Exam.integrity(score);
      score.increment(code.code(y, q), true);
      Exam.integrity(score);
      score.increment(code.code(z), false);
      Exam.integrity(score);
      score.increment(code.code(q), false);
      Exam.integrity(score);
    }

    final double expected = MathUtils.logFactorial(5) - MathUtils.LOG_2 * 5 - MathUtils.logFactorial(2) - MathUtils.logFactorial(3);
    //System.err.println(expected);
    assertEquals(expected, score.lnProbability(), 0.00000001);
  }
}
