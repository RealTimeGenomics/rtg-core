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

import com.rtg.util.intervals.Interval;
import com.rtg.util.intervals.RangeList;

import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

/**
 * Checks whether a record should be stripped according to given criteria. Provides stripping function.
 */
abstract class PositionAndStrandChecker {

  static final int MAX_OP_LEN = 20;

  protected final int mTolerance;
  protected final int[] mPosDiffStats;
  protected final int[] mSoftClipStats;
  protected final int[] mMismatchStats;
  protected final int[] mInsertStats;
  protected final int[] mDeletionStats;
  protected long mBasesTrimmed = 0;

  PositionAndStrandChecker(int tolerance) {
    mTolerance = tolerance;
    mPosDiffStats = new int[mTolerance * 2 + 1];
    mSoftClipStats = new int[MAX_OP_LEN];
    mMismatchStats = new int[MAX_OP_LEN];
    mInsertStats = new int[MAX_OP_LEN];
    mDeletionStats = new int[MAX_OP_LEN];
  }

  abstract boolean checkStrand(SAMRecord record);
  abstract boolean checkPosition(SAMRecord record, Interval data);

  abstract int getStartDataIndex(SAMRecord record, RangeList<?> list);

  abstract void stripRecord(SAMRecord record, SAMRecord mate, Interval data);

  protected void updateStrippedStats(CigarOperator operator, int consume) {
    final int statIndex = consume > MAX_OP_LEN ? MAX_OP_LEN - 1 : consume - 1;
    switch (operator) {
      case X:
        mMismatchStats[statIndex]++;
        break;
      case S:
        mSoftClipStats[statIndex]++;
        break;
      case I:
        mInsertStats[statIndex]++;
        break;
      case D:
        mDeletionStats[statIndex]++;
        break;
      default:
        break;
    }
  }

  static void updateTlenAndMateStart(SAMRecord rec1, SAMRecord rec2) {
    final SAMRecord leftMost;
    final SAMRecord rightMost;
    if (rec1.getAlignmentStart() <= rec2.getAlignmentStart()) {
      leftMost = rec1;
      rightMost = rec2;
    } else {
      leftMost = rec2;
      rightMost = rec1;
    }
    final int end = Math.max(leftMost.getAlignmentEnd(), rightMost.getAlignmentEnd());
    final int start = leftMost.getAlignmentStart() - 1;
    final int tlen = end - start;
    leftMost.setInferredInsertSize(tlen);
    rightMost.setInferredInsertSize(-tlen);
    leftMost.setMateAlignmentStart(rightMost.getAlignmentStart());
    rightMost.setMateAlignmentStart(leftMost.getAlignmentStart());
  }
}
