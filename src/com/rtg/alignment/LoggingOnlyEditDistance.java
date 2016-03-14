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
    mCount++;
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
