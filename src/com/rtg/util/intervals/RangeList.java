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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;

import com.rtg.util.diagnostic.NoTalkbackSlimException;

/**
 * A structure that holds multiple ranges that can be searched for a point within the ranges.
 *
 * Makes use of a segment tree <code>http://en.wikipedia.org/wiki/Segment_tree</code> to store and efficiently look up meta-data within ranges.
 * The input meta-data are split and merged to form a single set of non-overlapping meta-data ranges.
 * A continuous set of ranges (with and without meta-data) is stored in the ranges array.
 * A simple binary search mechanism is used to look up the meta-data for a given location.
 *
 */
public class RangeList<T> {

  /**
   * Range with meta-data class.
   */
  public static final class RangeData<T> extends Range {

    private List<T> mMeta = null;

    private List<RangeData<T>> mOriginalRanges = null;

    /**
     * Constructor.
     * @param range an existing range range.
     * @param meta the meta-data to store.
     */
    public RangeData(Interval range, T meta) {
      this(range.getStart(), range.getEnd(), meta);
    }

    /**
     * Constructor.
     * @param start the start of the range (0-based inclusive).
     * @param end the end of the range (0-based exclusive).
     * @param meta the meta-data to store.
     */
    public RangeData(int start, int end, T meta) {
      this(start, end, new ArrayList<T>());
      mMeta.add(meta);
    }

    RangeData(int start, int end) {
      this(start, end, (List<T>) null);
    }

    /**
     * Constructor.
     * @param start the start of the range (0-based inclusive).
     * @param end the end of the range (0-based exclusive).
     * @param metas the list of meta-data to store.
     */
    private RangeData(int start, int end, List<T> metas) {
      super(start, end);
      if (end < start) {
        throw new NoTalkbackSlimException("Invalid range : [" + start + "," + end + ") ");
      }
      mMeta = metas;
    }

    /**
     * Returns whether given location is within range
     *
     * @param loc sequence position
     * @return true if in range
     */
    public boolean isInRange(int loc) {
      return loc >= getStart() && loc < getEnd();
    }

    /**
     * Returns the meta information for this range
     *
     * @return meta list
     */
    public List<T> getMeta() {
      if (mMeta == null) {
        return null;
      }
      return new ArrayList<>(mMeta);
    }

    /**
     * Add single meta-data to range object.
     * @param meta the meta-data to add.
     */
    public void addMeta(T meta) {
      if (meta != null) {
        if (mMeta == null) {
          mMeta = new ArrayList<>();
        }
        mMeta.add(meta);
      }
    }

    /**
     * Add list of meta-data to range object.
     * @param metas the meta-data to add.
     */
    public void addMeta(List<T> metas) {
      if (metas != null) {
        if (mMeta == null) {
          mMeta = new ArrayList<>();
        }
        mMeta.addAll(metas);
      }
    }

    private void addOriginalRange(RangeData<T> originalRange) {
      if (mOriginalRanges == null) {
        mOriginalRanges = new ArrayList<>();
      }
      mOriginalRanges.add(originalRange);
    }

    /**
     * @return the list of ranges as originally specified that are covered by this range
     */
    public List<RangeData<T>> getOriginalRanges() {
      return mOriginalRanges;
    }

  }

  private final List<RangeData<T>> mRanges;
  private final List<RangeData<T>> mNonEmptyRanges;

  /**
   * Constructor.
   * @param ranges the list of meta-data ranges to store for searching.
   */
  public RangeList(List<RangeData<T>> ranges) {
    if (ranges == null || ranges.size() == 0) {
      mRanges = new ArrayList<>(1);
      mRanges.add(new RangeData<T>(Integer.MIN_VALUE, Integer.MAX_VALUE));
    } else {
      // get list of range boundaries
      final HashSet<Integer> pivots = new HashSet<>();
      for (final RangeData<T> range : ranges) {
        pivots.add(range.getStart());
        pivots.add(range.getEnd());
      }
      final int[] pivots2 = new int[pivots.size()];
      int i2 = 0;
      for (final Integer x : pivots) {
        pivots2[i2] = x;
        i2++;
      }
      Arrays.sort(pivots2);

      // set up continuous non-overlapping ranges for -inf to +inf
      mRanges = new ArrayList<>(pivots2.length + 1);
      if (pivots2[0] != Integer.MIN_VALUE) {
        mRanges.add(new RangeData<T>(Integer.MIN_VALUE, pivots2[0]));
      }
      for (int i = 1; i < pivots2.length; i++) {
        mRanges.add(new RangeData<T>(pivots2[i - 1], pivots2[i]));
      }
      if (pivots2[pivots2.length - 1] != Integer.MAX_VALUE) {
        mRanges.add(new RangeData<T>(pivots2[pivots2.length - 1], Integer.MAX_VALUE));
      }

      // load original meta into non-overlapping range structure
      for (final RangeData<T> range : ranges) {
        int index = findFullRangeIndex(range.getStart());
        while (index < mRanges.size() && mRanges.get(index).getEnd() <= range.getEnd()) {
          mRanges.get(index).addMeta(range.getMeta());
          mRanges.get(index).addOriginalRange(range);
          index++;
        }
      }
    }
    mNonEmptyRanges = new ArrayList<>();
    for (final RangeData<T> range : mRanges) {
      if (range.getMeta() != null) {
        mNonEmptyRanges.add(range);
      }
    }
  }

  /**
   * Return the list of meta-data for a given location.  Returns null if no meta-data exists.
   * @param loc the position to look up the containing range for.
   * @return meta-data list
   */
  public List<T> find(int loc) {
    return findRange(loc).getMeta();
  }

  private RangeData<T> findRange(int loc) {
    return mRanges.get(findFullRangeIndex(loc));
  }

  /**
   * @return the list of non-empty ranges
   */
  public List<RangeData<T>> getRangeList() {
    return mNonEmptyRanges;
  }

  /**
   * @return the full list of ranges, includes ranges corresponding to empty intervals
   */
  public List<RangeData<T>> getFullRangeList() {
    return mRanges;
  }

  /**
   * Return the index of the range within the full range list containing the specified point.
   * @param loc the position to search
   * @return the index of the range entry containing the position
   */
  public final int findFullRangeIndex(int loc) {
    if (loc == Integer.MAX_VALUE) {
      return mRanges.size() - 1;
    }
    int min = 0;
    int max = mRanges.size();
    int res = (max + min) / 2;

    // binary search
    boolean found = false;
    while (!found) {
      //System.err.println(min + " " + res + " " + max + " : " + loc + " " + ranges[res]);
      final RangeData<T> range = mRanges.get(res);
      if (range.isInRange(loc)) {
        found = true;
      } else {
        if (loc < range.getStart()) {
          max = res;
        } else {
          min = res;
        }
        res = (max + min) / 2;
      }
    }
    return res;
  }

  @Override
  public String toString() {
    return mNonEmptyRanges.toString();
  }
}
