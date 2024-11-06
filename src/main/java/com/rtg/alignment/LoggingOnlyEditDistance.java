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
package com.rtg.alignment;

import java.util.zip.CRC32;

import com.rtg.mode.DnaUtils;
import com.rtg.util.StringUtils;
import com.rtg.util.diagnostic.Diagnostic;

/**
 * An edit distance implementation that simply logs to the reads to the developer log.
 * Extremely useful in the <code>PrioritisedEditDistance</code> chain.
 *
 */
class LoggingOnlyEditDistance implements UnidirectionalEditDistance {

  private final int mMaxTimesToLog;
  private int mCount;
  private final String mLabel;


  /**
   * Creates a logging only edit distance. If the sub edit distance returns null then we log the read/template.
   * Only logs the <code>maxTimesToLog</code> calls.
   *
   * @param maxTimesToLog maximum reads to display
   */
  protected LoggingOnlyEditDistance(final int maxTimesToLog) {
    this(maxTimesToLog, "");
  }

  /**
   * Creates a logging only edit distance. If the sub edit distance returns null then we log the read/template.
   * Only logs the <code>maxTimesToLog</code> calls.
   *
   * @param maxTimesToLog maximum reads to display
   * @param label label to include in log
   */
  protected LoggingOnlyEditDistance(final int maxTimesToLog, final String label) {
    mMaxTimesToLog = maxTimesToLog;
    mCount = 0;
    mLabel = label;
  }

  /**
   * Creates a logging only edit distance. If the sub edit distance returns null then we log the read/template.
   * Only logs the first 3 calls.
   */
  LoggingOnlyEditDistance() {
    this(3, "");
  }

  @Override
  public int[] calculateEditDistance(byte[] read, int rlen, byte[] template, int zeroBasedStart, int maxScore, int maxShift, boolean cgLeft) {
    ++mCount;
    if (mCount <= mMaxTimesToLog) {
      final CRC32 checksum = new CRC32();
      checksum.update(read, 0, rlen);

      //      final GotohEditDistance g = new GotohEditDistance(1);
      //      final int[] gotoh = g.calculateEditDistance(read, rlen, template, zeroBasedStart, maxScore);
      //      sb.append("Gotoh AS: " + gotoh[ActionsHelper.ALIGNMENT_SCORE_INDEX]);

      Diagnostic.developerLog("LoggingOnlyEditDistance " + mLabel + " " + mCount + ": " + System.identityHashCode(read) + " crc=" + checksum.getValue() + " read.length=" + read.length + " rlen=" + rlen + StringUtils.LS + "zeroStart=" + zeroBasedStart + " maxScore=" + maxScore + " template.length=" + template.length + " cgLeft=" + cgLeft + StringUtils.LS + "ED: " + DnaUtils.bytesToSequenceIncCG(read) + " " + DnaUtils.bytesToSequenceIncCG(template, Math.max(0, zeroBasedStart - 20), rlen + 40) + " 20");
    }
    return null;
  }


  @Override
  public void logStats() {
    //empty
  }
  @Override
  public int[] calculateEditDistanceFixedBoth(byte[] read, int readStartPos, int readEndPos, byte[] template, int templateStartPos, int templateEndPos, int maxScore, int maxShift) {
    throw new UnsupportedOperationException();
  }
  @Override
  public int[] calculateEditDistanceFixedEnd(final byte[] read, final int readStartPos, final int readEndPos, final byte[] template, final int templateExpectedStartPos,
      final int templateEndPos, final int maxScore, int maxShift) {
    throw new UnsupportedOperationException();
  }
  @Override
  public int[] calculateEditDistanceFixedStart(byte[] read, int readStartPos, int readEndPos, byte[] template, int templateStartPos, int maxScore, int maxShift) {
    throw new UnsupportedOperationException();
  }
}
