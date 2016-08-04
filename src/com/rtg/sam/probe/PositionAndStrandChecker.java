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

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.util.intervals.RangeList;

import htsjdk.samtools.SAMRecord;

/**
 * Checks whether a record should be stripped according to given criteria. Provides stripping function.
 */
@TestClass("com.rtg.sam.probe.PosCheckerTest")
abstract class PositionAndStrandChecker {
  protected final int mTolerance;
  protected final int[] mStats;

  PositionAndStrandChecker(int tolerance) {
    mTolerance = tolerance;
    mStats = new int[mTolerance * 2 + 1];
  }

  abstract boolean check(SAMRecord record, RangeList.RangeData<String> data);

  abstract int getStartDataIndex(SAMRecord record, RangeList<String> list);

  abstract void stripRecord(SAMRecord record, RangeList.RangeData<String> data);
}
