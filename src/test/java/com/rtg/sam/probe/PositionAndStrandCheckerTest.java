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

package com.rtg.sam.probe;

import htsjdk.samtools.CigarOperator;
import junit.framework.TestCase;

public class PositionAndStrandCheckerTest extends TestCase {

  public void testStats() {
    final PosChecker pos = new PosChecker(10);
    assertEquals(10, pos.mTolerance);
    assertEquals(21, pos.mPosDiffStats.length);
    pos.updateStrippedStats(CigarOperator.D, 7);
    pos.updateStrippedStats(CigarOperator.D, 7);
    pos.updateStrippedStats(CigarOperator.I, 5);
    pos.updateStrippedStats(CigarOperator.S, 9);
    pos.updateStrippedStats(CigarOperator.X, 13);
    pos.updateStrippedStats(CigarOperator.D, 30);

    assertEquals(1, pos.mMismatchStats[12]);
    assertEquals(2, pos.mDeletionStats[6]);
    assertEquals(1, pos.mInsertStats[4]);
    assertEquals(1, pos.mSoftClipStats[8]);
    assertEquals(1, pos.mDeletionStats[19]);
  }
}
