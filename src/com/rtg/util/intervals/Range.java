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
 * A class representing a range of integers as a zero-based half-open interval.
 * May be used to represent a Locus, or other semantic range.
 */
public class Range implements Interval {

  private final int mStart;
  private final int mEnd;

  /**
   * @param start 0 based start position
   * @param end 0 based end position exclusive
   */
  public Range(int start, int end) {
    if (end < start) {
      throw new IllegalArgumentException("Locus start must be less than or equal to end. start=" + start + ", end=" + end);
      // When start == end, this denotes an insertion point between two bases
    }
    mStart = start;
    mEnd = end;
  }

  @Override
  public int getLength() {
    return getEnd() - getStart();
  }

  @Override
  public int getStart() {
    return mStart;
  }

  @Override
  public int getEnd() {
    return mEnd;
  }

  @Override
  public String toString() {
    return (getStart() + 1) + "-" + getEnd(); // Output as 1-based with end inclusive, like RegionRestriction
  }

}
