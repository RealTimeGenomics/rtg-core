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
import com.rtg.variant.bayes.Code;
import com.rtg.variant.bayes.CodeDiploid;

import junit.framework.TestCase;

/**
 */
public class SegregationScoreTest extends TestCase {

  public void test() {
    final Code code = new CodeDiploid(2);
    final SegregationScore score = new SegregationScore(code, code.code(0, 1), code.code(0, 1));
    assertTrue(score.integrity());
    assertEquals(0.0, score.lnProbability());

    score.increment(code.code(0, 0), true);
    assertTrue(score.integrity());
    score.increment(code.code(0, 0), true);
    assertTrue(score.integrity());
    score.increment(code.code(0, 0), true);
    assertTrue(score.integrity());

    score.increment(code.code(0, 1), true);
    assertTrue(score.integrity());
    score.increment(code.code(1, 0), true);
    assertTrue(score.integrity());
    final double expected = Math.log(120) + Math.log(0.5) * 8 - MathUtils.logFactorial(3) - MathUtils.logFactorial(2);
    //System.err.println(expected);
    assertEquals(expected, score.lnProbability(), 0.000001);
  }

  public void test2() {
    final Code code = new CodeDiploid(3);
    final SegregationScore score = new SegregationScore(code, code.code(0, 1), code.code(0, 1));
    assertTrue(score.integrity());
    assertEquals(0.0, score.lnProbability());

    score.increment(code.code(1, 2), true); // Test having a code not in parents
    assertTrue(score.integrity());

    score.increment(code.code(0, 0), true);
    assertTrue(score.integrity());
    score.increment(code.code(0, 1), true);
    assertTrue(score.integrity());
    score.increment(code.code(1, 1), true);
    assertTrue(score.integrity());

    final double expected = Math.log(6) + Math.log(0.5) * 5;
    //System.err.println(expected);
    assertEquals(expected, score.lnProbability(), 0.000001);
  }

}
