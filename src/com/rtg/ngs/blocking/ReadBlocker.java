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
package com.rtg.ngs.blocking;

import java.io.Closeable;

import com.rtg.util.License;
import com.rtg.util.diagnostic.Diagnostic;

/**
 * Maintain a count of records for identifiers.
 *
 */
public class ReadBlocker implements Closeable {

  private static final int MAX_COUNT = Short.MAX_VALUE * 2 + 1;
  private final short[] mCounts; // treated as unsigned here
  private final int mThreshold;
  private final String mTitle;

  /**
  * Copy constructor for making a non-synchronized version from a synchronized one
   * @param source the source ReadBlocker
   */
  public ReadBlocker(ReadBlocker source) {
    mTitle = source.mTitle;
    mCounts = source.mCounts;
    mThreshold = source.mThreshold;
  }

  /**
   * Creates a counter for <code>count</code> records blocking at <code>
   * threshold</code>.
   *
   * @param count number of reads
   * @param threshold blocking threshold in range 1 to 255
   * @param title a title to use during logging
   */
  public ReadBlocker(final long count, final int threshold, final String title) {
    if (count > Integer.MAX_VALUE) {
      throw new IllegalArgumentException("Too many reads");
    }
    if (threshold != -1 && (threshold > MAX_COUNT || threshold < 1)) {
      throw new IllegalArgumentException();
    }
    mCounts = new short[(int) count];
    mThreshold = threshold;
    mTitle = title;
  }

  /**
   * Creates a counter for <code>count</code> records blocking at <code>
   * threshold</code>.
   *
   * @param count number of reads
   * @param threshold blocking threshold in range 1 to 255.
   */
  public ReadBlocker(final long count, final int threshold) {
    this(count, threshold, "blocked pairings");
  }

  /**
   * Resets the count for the given read.
   *
   * @param r read number
   */
  public void reset(final int r) {
    mCounts[r] = 0;
  }

  /**
   * Increment the count for the given read.
   *
   * @param r read number
   */
  public void increment(final int r) {
    if (mCounts[r] != -1) {
      mCounts[r]++;
    }
  }

  /**
   * Check if the specified read is blocked.
   *
   * @param r read to check
   * @return true if blocked
   */
  public boolean isBlocked(final int r) {
    if (mThreshold == -1) {
      return false;
    }
    return (mCounts[r] & MAX_COUNT) >= mThreshold;
  }

  @Override
  public void close() {
    if (License.isDeveloper()) {
      final int[] h = new int[MAX_COUNT + 1];
      for (final int k : mCounts) {
        h[k & MAX_COUNT]++;
      }
      Diagnostic.developerLog("Statistics of " + mTitle);
      long sum = 0;
      for (int k = 0; k < MAX_COUNT + 1; ++k) {
        if (h[k] > 0) {
          sum += h[k];
          final String c = k == MAX_COUNT ? ">= " + MAX_COUNT : String.valueOf(k);
          Diagnostic.developerLog(h[k] + " reads had count " + c);
        }
      }
      Diagnostic.developerLog("Total reads " + sum);
    }
  }

  /**
   * The number of records with this read identifier that have been written.
   * If this is equal to or greater than the threshold, then it means
   * that a large/infinite number of records have been written.
   * @param r Read identifier
   * @return the number of records written, or the threshold.
   */
  public int getCount(int r) {
    return mCounts[r] & MAX_COUNT;
  }
}
