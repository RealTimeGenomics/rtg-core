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

import java.util.Map;
import java.util.TreeMap;

import com.reeltwo.jumble.annotations.TestClass;

/**
 * Store information about multiple intervals on a single sequence.
 * Overlapping regions are joined together into a single interval.
 */
@TestClass("com.rtg.util.intervals.ReferenceRegionsTest")
class MergedIntervals {

  /**
   * Map of region start to region end
   * a position is within the region if the floor entries value is greater than it
   */
  final TreeMap<Integer, Integer> mIntervals = new TreeMap<>();

  /**
   * Add a new interval to the set
   * @param interval the region to add
   */
  void add(Interval interval) {
    add(interval.getStart(), interval.getEnd());
  }

  /**
   * Add a new region to the set
   * @param start 0-based inclusive start position of the region
   * @param end 0-based exclusive end position of the region
   */
  void add(int start, int end) {
    final Map.Entry<Integer, Integer> floor = mIntervals.floorEntry(start);
    final Map.Entry<Integer, Integer> endFloor = mIntervals.floorEntry(end);
    final int actualStart;
    final int actualEnd;
    if (floor != null && start >= floor.getKey() && end <= floor.getValue()) {
      //Don't bother with regions already fully contained
      return;
    }
    //Get Start for Region
    if (floor == null || start > floor.getValue()) {
      actualStart = start;
    } else {
      actualStart = Math.min(floor.getKey(), start);
    }
    //Get End for Region
    if (endFloor == null) {
      actualEnd = end;
    } else {
      actualEnd = Math.max(endFloor.getValue(), end);
    }
    //Remove any existing regions contained within the new bounds
    Map.Entry<Integer, Integer> overlappedRegion;
    while ((overlappedRegion = mIntervals.floorEntry(actualEnd)) != null) {
      if (overlappedRegion.getKey() < actualStart) {
        break;
      }
      mIntervals.remove(overlappedRegion.getKey());
    }
    mIntervals.put(actualStart, actualEnd);
  }

  /**
   * @param pos zero based position within the sequence
   * @return true if the position provided falls within the intervals
   */
  boolean enclosed(int pos) {
    final Map.Entry<Integer, Integer> floor = mIntervals.floorEntry(pos);
    return floor != null && floor.getValue() > pos;
  }

  /**
   * @param start zero based position within the sequence
   * @param end zero based position within the sequence
   * @return true if the position provided falls entirely within the intervals
   */
  boolean enclosed(int start, int end) {
    final Map.Entry<Integer, Integer> floor = mIntervals.floorEntry(end);
    if (floor == null) {
      return false;
    } else if (end > floor.getValue()) {
      return false;
    } else if (start < floor.getKey()) {
      return false;
    }
    return true;
  }

  /**
   * @param start zero based start position within the sequence
   * @param end zero based end position within the sequence
   * @return true if the range specified is overlapped by the intervals
   */
  boolean overlapped(int start, int end) {
    final Map.Entry<Integer, Integer> floor = mIntervals.floorEntry(end);
    if (floor == null) {
      return false;
    } else {
      return start < floor.getValue() && end > floor.getKey();
    }
  }

  /**
   * Work out the total length of all intervals
   * @return the total length covered by regions
   */
  public int totalLength() {
    int total = 0;
    for (Map.Entry<Integer, Integer> region : mIntervals.entrySet()) {
      total += region.getValue() - region.getKey();
    }
    return total;
  }
}
