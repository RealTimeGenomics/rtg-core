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
