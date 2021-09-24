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
