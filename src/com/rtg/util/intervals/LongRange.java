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

/**
 * A class representing a range of non-negative longs as a zero-based half-open interval.
 */
public class LongRange {

  /** Indicates either the start or end is missing (usually meaning the range should be extended as far as possible) */
  public static final long MISSING = -1; // This should probably actually be MIN or MAX, as it interferes with the ability to use LongRange for negative intervals, but currently too much depends on hardcoded -1s.

  /** A singleton that treats all locations as in range. */
  public static final LongRange NONE = new LongRange(MISSING, MISSING) {  // Long.MIN_VALUE, Long.MAX_VALUE) {
    @Override
    public boolean isInRange(final long value) {
      return true;
    }
    @Override
    public String toString() {
      return "[All inclusive]";
    }
  };

  private final long mStart;
  private final long mEnd;

  /**
   * @param start 0 based start position
   * @param end 0 based end position exclusive
   */
  public LongRange(long start, long end) {
    if (end != MISSING && start != MISSING && end < start) {
      throw new IllegalArgumentException("Locus start must be less than or equal to end. start=" + start + ", end=" + end);
      // When start == end, this denotes an insertion point between two bases
    }
    mStart = start;
    mEnd = end;
  }

  /**
   * Nonnegative length of this interval.
   *
   * @return length
   */
  public long getLength() {
    return getEnd() - getStart();
  }

  /**
   * Zero-based start of the range.
   * @return start of the range
   */
  public long getStart() {
    return mStart;
  }

  /**
   * Zero-based end of the range (exclusive).
   * @return end of the range
   */
  public long getEnd() {
    return mEnd;
  }

  /**
   * @param value the query position
   * @return true if the query position is within the range region.
   */
  public boolean isInRange(final long value) {
    return value >= mStart && value < mEnd;
  }

  @Override
  public String toString() {
    return (getStart() + 1) + "-" + getEnd(); // Output as 1-based with end inclusive, like RegionRestriction
  }

}
