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
    first.mReadThroughOnR1 = 10;
    first.mReadThroughOnR2 = 20;
    first.mMajorOverlap = 30;
    first.mPartialOverlap = 40;
    first.mNoAlignment = 50;
    first.mPoorAlignment = 60;
    first.mTotal = 70;
    second.mReadThroughOnR1 = 100;
    second.mReadThroughOnR2 = 200;
    second.mMajorOverlap = 300;
    second.mPartialOverlap = 400;
    second.mNoAlignment = 500;
    second.mPoorAlignment = 600;
    second.mTotal = 700;

    accumulator.accumulate(first);
    accumulator.accumulate(second);

    assertEquals(110, accumulator.mReadThroughOnR1);
    assertEquals(220, accumulator.mReadThroughOnR2);
    assertEquals(330, accumulator.mMajorOverlap);
    assertEquals(440, accumulator.mPartialOverlap);
    assertEquals(550, accumulator.mNoAlignment);
    assertEquals(660, accumulator.mPoorAlignment);
    assertEquals(770, accumulator.mTotal);

  }

}