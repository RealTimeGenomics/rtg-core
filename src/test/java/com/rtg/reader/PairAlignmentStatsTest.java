/*
 * Copyright (c) 2017. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
 */

package com.rtg.reader;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

/**
 */
public class PairAlignmentStatsTest {
  @Test
  public void accumulateSumsAllFields() {
    final PairAlignmentStats accumulator = new PairAlignmentStats();
    final PairAlignmentStats first = new PairAlignmentStats();
    final PairAlignmentStats second = new PairAlignmentStats();
    first.mR1ReadThrough = 10;
    first.mR2ReadThrough = 20;
    first.mR2ReadIntoR1Probe = 30;
    first.mR1ReadIntoR2Probe = 10;
    first.mOverlapping = 40;
    first.mNoAlignment = 50;
    first.mPoorAlignment = 60;
    first.mTotalInput = 70;
    first.mTotalOutput = 69;
    second.mR1ReadThrough = 100;
    second.mR2ReadThrough = 200;
    second.mR2ReadIntoR1Probe = 300;
    second.mR1ReadIntoR2Probe = 100;
    second.mOverlapping = 400;
    second.mNoAlignment = 500;
    second.mPoorAlignment = 600;
    second.mTotalInput = 700;
    second.mTotalOutput = 695;

    accumulator.accumulate(first);
    accumulator.accumulate(second);

    assertEquals(110, accumulator.mR1ReadThrough);
    assertEquals(220, accumulator.mR2ReadThrough);
    assertEquals(330, accumulator.mR2ReadIntoR1Probe);
    assertEquals(110, accumulator.mR1ReadIntoR2Probe);
    assertEquals(440, accumulator.mOverlapping);
    assertEquals(550, accumulator.mNoAlignment);
    assertEquals(660, accumulator.mPoorAlignment);
    assertEquals(770, accumulator.mTotalInput);
    assertEquals(764, accumulator.mTotalOutput);

  }

}