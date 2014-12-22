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
package com.rtg.tabix;

import java.util.Arrays;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.util.QuickSort;
import com.rtg.util.intervals.SequenceNameLocus;

/**
 * Hold pairs of start and end virtual file offsets as looked up from a tabix or BAM index. Each pair may have a region associated with it.
 */
@TestClass("com.rtg.tabix.TabixIndexReaderTest")
public class VirtualOffsets {
  private static final long[] EMPTY_POS = {};
  private static final SequenceNameLocus[] EMPTY_LOCUS = {};

  private long[] mStarts;
  private long[] mEnds;
  private SequenceNameLocus[] mRegions;
  private int mSize;

  /**
   * Convenience constructor for an initially empty offsets list.
   */
  public VirtualOffsets() {
    this(EMPTY_POS, EMPTY_POS, EMPTY_LOCUS);
  }

  /**
   * Convenience constructor for a single offset pair.
   * @param start start virtual offset
   * @param end end virtual offset
   * @param region the region corresponding to the offset pair.
   */
  public VirtualOffsets(long start, long end, SequenceNameLocus region) {
    this(new long[]{start}, new long[] {end}, new SequenceNameLocus[] {region});
  }

  /**
   * Create the offset storage.
   * @param starts the set of start virtual offsets
   * @param ends the set of end virtual offsets
   * @param regions the set of regions corresponding to each start/end pair
   */
  public VirtualOffsets(long[] starts, long[] ends, SequenceNameLocus[] regions) {
    if (starts == null || ends == null || regions == null) {
      throw new NullPointerException();
    }
    if (starts.length != ends.length) {
      throw new IllegalArgumentException("Start and end virtual offset lists must have same length");
    }
    if (starts.length != regions.length) {
      throw new IllegalArgumentException("Position and region lists must have same length");
    }
    mStarts = starts;
    mEnds = ends;
    mRegions = regions;
    mSize = mStarts.length;
  }

  /**
   * @return the number of virtual offset pairs
   */
  public int size() {
    return mSize;
  }

  /**
   * Adds a new start, end pair.
   * @param start the start virtual offset
   * @param end the end virtual offset
   * @param region the region corresponding to the offset pair.
   */
  public void add(long start, long end, SequenceNameLocus region) {
    if (mSize == mStarts.length) {
      final int newLen = mStarts.length * 3 / 2 + 1;
      mStarts = Arrays.copyOf(mStarts, newLen);
      mEnds = Arrays.copyOf(mEnds, newLen);
      mRegions = Arrays.copyOf(mRegions, newLen);
    }
    mStarts[mSize] = start;
    mEnds[mSize] = end;
    mRegions[mSize] = region;
    mSize++;
  }

  /**
   * Gets the start offset of the specified pair
   * @param index the pair index
   * @return the start virtual offset
   */
  public long start(int index) {
    return mStarts[index];
  }

  /**
   * Gets the end offset of the specified pair
   * @param index the pair index
   * @return the end virtual offset
   */
  public long end(int index) {
    return mEnds[index];
  }

  /**
   * Gets the region of the specified pair
   * @param index the pair index
   * @return the region
   */
  public SequenceNameLocus region(int index) {
    return mRegions[index];
  }

  /**
   * Sort all offsets by increasing starting position.
   */
  public void sort() {
    QuickSort.sort(new QuickSortProxy());
  }

  /**
   * @param l raw virtual offset
   * @return slightly more human readable format of virtual offset
   */
  public static String offsetToString(long l) {
    return "(" + (l >>> 16) + ", " + (l & 0xFFFF) + ")";
  }

  @Override
  public String toString() {
    final StringBuilder sb = new StringBuilder("[");
    for (int i = 0; i < size(); i++) {
      if (i > 0) {
        sb.append(", ");
      }
      sb.append("{").append(offsetToString(start(i))).append(", ").append(offsetToString(end(i)));
      if (region(i) != null) {
        sb.append(", ").append(region(i));
      }
      sb.append("}");
    }
    sb.append("]");
    return sb.toString();
  }


  /**
   * Quick sort proxy
   */
  public class QuickSortProxy implements QuickSort.SortProxy {

    @Override
    public int compare(long index1, long index2) {
      final int res = Long.compare(mStarts[(int) index1], mStarts[(int) index2]);
      if (res != 0) {
        return res;
      }
      final SequenceNameLocus region1 = mRegions[(int) index1];
      final SequenceNameLocus region2 = mRegions[(int) index2];
      if (region1 == null || region2 == null) { // Only if one of the regions corresponds to no restriction
        return 0;
      }
      // This sorting assumes non-overlapping regions
      if (region1.getSequenceName().equals(region2.getSequenceName())) {
        // Two regions in the same block and on the same sequence, want the lowest coord first
        return Integer.compare(region1.getStart(), region2.getStart());
      } else {
        // Two regions in the same block and but on different sequences, want the highest coord (end of sequence) before lowest coord (start of next)
        return Integer.compare(region2.getStart(), region1.getStart());
      }
    }

    @Override
    public long length() {
      return mSize;
    }

    @Override
    public void swap(long index1, long index2) {
      final long t = mStarts[(int) index1];
      mStarts[(int) index1] = mStarts[(int) index2];
      mStarts[(int) index2] = t;

      final long t2 = mEnds[(int) index1];
      mEnds[(int) index1] = mEnds[(int) index2];
      mEnds[(int) index2] = t2;

      final SequenceNameLocus t3 = mRegions[(int) index1];
      mRegions[(int) index1] = mRegions[(int) index2];
      mRegions[(int) index2] = t3;
    }

  }
}
