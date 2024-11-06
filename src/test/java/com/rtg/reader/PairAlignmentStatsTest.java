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
