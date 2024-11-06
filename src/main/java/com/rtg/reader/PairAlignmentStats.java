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
