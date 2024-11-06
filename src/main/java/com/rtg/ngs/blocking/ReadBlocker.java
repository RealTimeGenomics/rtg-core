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
