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

import com.rtg.util.Histogram;
import com.rtg.util.TextTable;
import com.rtg.util.Utils;
import com.rtg.util.diagnostic.Diagnostic;

/**
 * Pair-alignment stats
 */
class PairAlignmentStats {
  protected long mR2ReadThrough;
  protected long mR1ReadThrough;
  protected long mR2ReadIntoR1Probe;
  protected long mR1ReadIntoR2Probe;
  protected long mOverlapping;
  protected long mNoAlignment;
  protected long mPoorAlignment;
  protected long mTotalInput;
  protected long mTotalOutput;
  protected Histogram mFragLengths = new Histogram();
  protected Histogram mOverlapDist = new Histogram();

  protected void accumulate(PairAlignmentStats other) {
    mR1ReadThrough += other.mR1ReadThrough;
    mR2ReadThrough += other.mR2ReadThrough;
    mR2ReadIntoR1Probe += other.mR2ReadIntoR1Probe;
    mR1ReadIntoR2Probe += other.mR1ReadIntoR2Probe;
    mOverlapping += other.mOverlapping;
    mNoAlignment += other.mNoAlignment;
    mPoorAlignment += other.mPoorAlignment;
    mTotalInput += other.mTotalInput;
    mTotalOutput += other.mTotalOutput;
    mFragLengths.addHistogram(other.mFragLengths);
    mOverlapDist.addHistogram(other.mOverlapDist);
    Diagnostic.userLog(toString());
  }

  protected String printSummary() {
    final TextTable t = new TextTable();
    t.addRow("Total input pairs", String.valueOf(mTotalInput), "");
    t.addRow("No alignment", String.valueOf(mNoAlignment), perc(mNoAlignment, mTotalInput));
    t.addRow("Poor alignment", String.valueOf(mPoorAlignment), perc(mPoorAlignment, mTotalInput));
    t.addRow("Overlapping", String.valueOf(mOverlapping), perc(mOverlapping, mTotalInput));
    t.addRow("R1 read through", String.valueOf(mR1ReadThrough), perc(mR1ReadThrough, mTotalInput));
    t.addRow("R2 read through", String.valueOf(mR2ReadThrough), perc(mR2ReadThrough, mTotalInput));
    if (mR1ReadIntoR2Probe > 0 || mR2ReadIntoR1Probe > 0) {
      t.addRow("R1 read into R2 probe", String.valueOf(mR1ReadIntoR2Probe), perc(mR1ReadIntoR2Probe, mTotalInput));
      t.addRow("R2 read into R1 probe", String.valueOf(mR2ReadIntoR1Probe), perc(mR2ReadIntoR1Probe, mTotalInput));
    }
    t.addRow("Total output pairs", String.valueOf(mTotalOutput), perc(mTotalOutput, mTotalInput));
    return t.toString();
  }

  @Override
  public String toString() {
    return "Reads-in: " + mTotalInput
      + " No-alignment: " + mNoAlignment
      + " Poor-alignment: " + mPoorAlignment
      + " Overlapping: " + mOverlapping
      + " R1-read-through: " + mR1ReadThrough
      + " R2-read-through: " + mR2ReadThrough
      + " R1-read-into-R2-probe: " + mR1ReadIntoR2Probe
      + " R2-read-into-R1-probe: " + mR2ReadIntoR1Probe
      + " Reads-out: " + mTotalOutput;
  }

  private String perc(long num, long total) {
    return Utils.realFormat(100.0 * num / total, 2) + "%";
  }
}
