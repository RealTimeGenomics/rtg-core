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
package com.rtg.index;

import com.rtg.util.array.CommonIndex;

/**
 */
public final class SearchUtils {

  /** Prevent instantiation. */
  private SearchUtils() {
  }

  /**
   * Searches part of a <code>CommonIndex</code> for the specified value.
   * Performs a binary search, so the array must be sorted.
   *
   * @param a the array to search
   * @param lowParam minimum index to search in
   * @param highParam maximum index to search in
   * @param key value to look for
   * @return The index of <code>key</code> if it is found, otherwise
   *         <code>-(<i>insertion point</i> + 1)</code>. (If found, the
   *         return value will always be <code>&gt;= 0</code>.)
   */
  public static long binarySearch(final CommonIndex a, final long lowParam, final long highParam, final long key) {
    long low = lowParam;
    long high = highParam;
    while (low <= high) {
      final long mid = (low + high) / 2;
      final long midVal = a.get(mid);
      if (midVal < key) {
        low = mid + 1;
      } else if (midVal > key) {
        high = mid - 1;
      } else {
        return mid; // key found
      }
    }
    return -(low + 1); // key not found.
  }

  /**
   * Searches part of a <code>CommonIndex</code> for the specified value.
   * Performs a search which alternates interpolation and binary searches, so the array must be sorted.
   *
   * @param a the array to search
   * @param lowParam minimum index to search in
   * @param highParam maximum index to search in
   * @param key value to look for
   * @return The index of <code>key</code> if it is found, otherwise
   *         <code>-(<i>insertion point</i> + 1)</code>. (If found, the
   *         return value will always be <code>&gt;= 0</code>.)
   */
  public static long interpolationSearch(final CommonIndex a, final long lowParam, final long highParam, final long key) {
    //System.err.println("key=" + key);
    long low = lowParam;
    long high = highParam;
    long lowf = a.get(low);
    if (key < lowf) {
      return -(low + 1);
    }
    long highf = a.get(high);
    if (key > highf) {
      return -(high + 1);
    }
    boolean interpolate = true;
    while (low <= high) {
      //System.err.println("low=" + low + " high=" + high + " lowf=" + lowf + " highf=" + highf);
      assert lowf <= key && key <= highf; // : lowf + ":" + key + ":" + highf;
      final long mid;
      if (interpolate) {
        final double frac = (key - lowf) / (double) (highf - lowf);
        final long length = high - low;
        mid = (long) (frac * length) + low;
      } else {
        mid = (low + high) / 2;
      }
      final long midVal = a.get(mid);
      //System.err.println("mid=" + mid + " midVal=" + midVal);
      if (midVal < key) {
        low = mid + 1;
        lowf = midVal + 1;
      } else if (midVal > key) {
        high = mid - 1;
        highf = midVal - 1;
      } else {
        return mid; // key found
      }
      interpolate = !interpolate;
    }
    return -(low + 1); // key not found.
  }

  /**
   * Searches part of a <code>CommonIndex</code> for adjacent values that bracket the key.
   * The returned value r satisfies <code>a[r] &le; key &lt; a[r+1]</code> or if no bracketing pair can be found
   * then returns -1.
   * Performs a search which alternates interpolation and binary searches, so the array must be sorted.
   *
   * @param a the array to search
   * @param lowParam minimum index to search in
   * @param highParam maximum index to search in
   * @param key value to look for
   * @return The low position of a pair of adjacent bracketing entries.
   */
  public static long bracketSearch(final CommonIndex a, final long lowParam, final long highParam, final long key) {
    //System.err.println("key=" + key);
    long low = lowParam;
    long high = highParam;
    long lowf = a.get(low);
    if (key < lowf) {
      return -1;
    }
    long highf = a.get(high);
    if (key >= highf) {
      return -1;
    }
    boolean interpolate = false;
    while (true) {
      if (low + 1 == high) {
        return low;
      }
      if (low >= high) { //key not found
        return -1;
      }
      //System.err.println("low=" + low + " high=" + high + " lowf=" + lowf + " highf=" + highf);
      assert lowf <= key && key < highf; // : lowf + ":" + key + ":" + highf;
      final long mid;
      if (interpolate) {
        final double frac = (key - lowf) / (double) (highf - lowf);
        final double frac1 = 0.9 * (frac - 0.5) + 0.5;
        final long length = high - low;
        mid = (long) (frac1 * length) + low;
      } else {
        mid = (low + high) / 2;
      }
      final long midVal = a.get(mid);
      //System.err.println("mid=" + mid + " midVal=" + midVal);
      if (midVal <= key) {
        low = mid;
        lowf = midVal;
      } else {
        high = mid;
        highf = midVal;
      }
      interpolate = !interpolate;
    }
  }


}
