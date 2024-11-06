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
package com.rtg.position.output;

import junit.framework.TestCase;

/**
 * Tests for Long mapping Mismatch scoring
 */
public class MismatchScoresTest extends TestCase {


  public void test() {

    final int maxGap = 40;
    final int maxIndel = 7;
    final int wordSize = 12;
    final int stepSize = 5;
    final MismatchScores ms = new MismatchScores(maxGap, maxIndel, wordSize, stepSize);
    ms.globalIntegrity();

    assertEquals(3, ms.minSubs(25));
    assertEquals(3, ms.minSubs(26));

    assertEquals(-7, ms.minDelta());
    assertEquals(+7, ms.maxDelta());

    assertEquals(-3.0, ms.score(0, 26, 0, 26));
    assertEquals(-4.0, ms.score(0, 26, 0, 27));

    assertEquals(2, ms.minSubs(20));
    assertEquals(-8.0, ms.score(0, 21, 0, 17));


    assertEquals(Double.NEGATIVE_INFINITY, ms.score(0, 30, 0, 40));
    assertEquals(Double.NEGATIVE_INFINITY, ms.score(0, 40, 0, 41));
    assertEquals("MismatchProbabilities", ms.toString());
  }

  public void testScoreMax() {
    final GapScorer ms = new MismatchScores(40, 7, 12, 5);
    assertEquals(ms.score(6, -1, 6, -1), ms.scoreMax(6, -1, -1, -1));
  }

  public void testMinSubs1() {
    final int maxGap = 20;
    final int maxIndel = 0;
    final int wordSize = 3;
    final int stepSize = 2;
    final MismatchScores ms = new MismatchScores(maxGap, maxIndel, wordSize, stepSize);
    ms.globalIntegrity();
    assertEquals(2, ms.minSubs(7));
  }

  public void testMinSubs2() {
    final int maxGap = 20;
    final int maxIndel = 0;
    final int wordSize = 5;
    final int stepSize = 2;
    final MismatchScores ms = new MismatchScores(maxGap, maxIndel, wordSize, stepSize);
    ms.globalIntegrity();
    assertEquals(2, ms.minSubs(7));
  }

  public void testMinSubs3() {
    final int maxGap = 1;
    final int maxIndel = 0;
    final int wordSize = 4;
    final int stepSize = 4;
    final MismatchScores ms = new MismatchScores(maxGap, maxIndel, wordSize, stepSize);
    ms.globalIntegrity();
    assertEquals(0, (int) ms.score(0, -1, 0, -1));
  }
}
