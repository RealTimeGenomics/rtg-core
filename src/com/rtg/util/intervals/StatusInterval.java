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
package com.rtg.util.intervals;

import java.util.Arrays;

/**
 * Finite width interval containing status values. This is backed by an array so memory issue may be a concern for large intervals.
 */
public class StatusInterval {

  /** The byte used to represent absence from the interval */
  public static final byte EMPTY = 0; // If you change this, the c'tor will need an explicit fill.

  private final byte[] mInterval;
  private final int mStart;

  /**
   * Construct a new interval spanning a given region.
   * @param start left edge of interval
   * @param end exclusive right edge of interval
   */
  public StatusInterval(final int start, final int end) {
    if (end <= start) {
      throw new IllegalArgumentException();
    }
    mStart = start;
    mInterval = new byte[end - start];
  }

  /**
   * Add points to the interval. Regions outside the overall interval are ignored.
   * @param start start position
   * @param end end position, exclusive
   * @param status the status to associate with points in the interval.
   */
  public void add(final int start, final int end, final byte status) {
    final int s = start - mStart;
    if (s < mInterval.length && end > mStart) {
      Arrays.fill(mInterval, Math.max(s, 0), Math.min(end - mStart, mInterval.length), status);
    }
  }

  /**
   * Test if given point is has a non-empty status in this interval.
   * @param pos position to test
   * @return true if point is in the interval
   */
  public boolean contains(final int pos) {
    final int index = pos - mStart;
    return mInterval[index] != EMPTY;
  }

  /**
   * Get the status value at a given point is contained in this interval.
   * @param pos position to test
   * @return the value stored in the interval. EMPTY is used for positions not present.
   */
  public byte get(final int pos) {
    final int index = pos - mStart;
    return mInterval[index];
  }
}
