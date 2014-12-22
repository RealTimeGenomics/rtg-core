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
