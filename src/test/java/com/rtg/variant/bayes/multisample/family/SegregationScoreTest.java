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
