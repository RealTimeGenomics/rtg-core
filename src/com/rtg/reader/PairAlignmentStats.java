/*
 * Copyright (c) 2016. Real Time Genomics Limited.
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

import com.rtg.util.Utils;
import com.rtg.util.diagnostic.Diagnostic;

/**
 * Pair-alignment stats
 */
class PairAlignmentStats {
  protected int mReadThroughOnR2;
  protected int mReadThroughOnR1;
  protected int mMajorOverlap;
  protected int mPartialOverlap;
  protected int mNoAlignment;
  protected int mPoorAlignment;
  protected int mTotal;

  protected void accumulate(PairAlignmentStats other) {
    mReadThroughOnR1 += other.mReadThroughOnR1;
    mReadThroughOnR2 += other.mReadThroughOnR2;
    mMajorOverlap += other.mMajorOverlap;
    mPartialOverlap += other.mPartialOverlap;
    mNoAlignment += other.mNoAlignment;
    mPoorAlignment += other.mPoorAlignment;
    mTotal += other.mTotal;
    printSummary();
  }

  protected void printSummary() {
    Diagnostic.userLog("Reads: " + mTotal
      + " No-alignment: " + perc(mNoAlignment, mTotal)
      + " Poor-alignment: " + perc(mPoorAlignment, mTotal)
      + " Partial-overlap: " + perc(mPartialOverlap, mTotal)
      + " Major-overlap: " + perc(mMajorOverlap, mTotal)
      + " Read-through-on-R1: " + perc(mReadThroughOnR1, mTotal)
      + " Read-through-on-R2: " + perc(mReadThroughOnR2, mTotal));
  }

  private String perc(int num, int total) {
    return Utils.realFormat(100.0 * num / total, 3);
  }
}
