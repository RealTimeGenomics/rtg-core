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
