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

package com.rtg.util;

import java.util.List;

/**
 * Attempt at making comparator methods easier to write.
 * Chain calls to compare and take the result at the end.
 */
public final class CompareHelper {
  private int mCurrent = 0;

  /**
   * Call this once for each pair of values you would like to sort by
   * Result will be the first non-zero comparison
   * @param o1 first comparable object
   * @param o2 object to compare to
   * @param <T> the type of objects to compare
   * @return this for chaining purposes
   */
  public <T extends Comparable<T>> CompareHelper compare(T o1, T o2) {
    if (mCurrent != 0) {
      return this;
    }
    mCurrent = o1.compareTo(o2);
    return this;
  }

  /**
   * Sets the comparison result to the first non-zero comparison of elements at the same index in the lists
   * or compares the sizes if one list is a prefix of the other
   * @param list1 first set of elements to compare
   * @param list2 second set of elements to compare
   * @param <T> type of the list
   * @return this for chaining purposes
   */
  public <T extends Comparable<T>> CompareHelper compareList(List<T> list1, List<T> list2) {
    if (mCurrent != 0) {
      return this;
    }
    for (int i = 0; i < list1.size() && i < list2.size(); ++i) {
      mCurrent = list1.get(i).compareTo(list2.get(i));
      if (mCurrent != 0) {
        return this;
      }
    }
    mCurrent = Integer.compare(list1.size(), list2.size());
    return this;
  }

  /**
   * Sets the comparison result to the first non-zero comparison of elements at the same index in the lists
   * or compares the sizes if one list is a prefix of the other
   * @param list1 first set of elements to compare
   * @param list2 second set of elements to compare
   * @param <T> type of the list
   * @return this for chaining purposes
   */
  public <T extends Comparable<T>> CompareHelper compareArray(T[] list1, T[] list2) {
    if (mCurrent != 0) {
      return this;
    }
    for (int i = 0; i < list1.length && i < list2.length; ++i) {
      mCurrent = list1[i].compareTo(list2[i]);
      if (mCurrent != 0) {
        return this;
      }
    }
    mCurrent = Integer.compare(list1.length, list2.length);
    return this;
  }

  /**
   * Add a result of an externally called comparison.
   * @param external result of some external attribute comparison
   * @return this for chaining purposes
   */
  public CompareHelper compare(int external) {
    if (mCurrent != 0) {
      return this;
    }
    mCurrent = external;
    return this;
  }

  /**
   * @return the value of the first comparison that was not zero or zero if all comparisons were equal
   */
  public int result() {
    return mCurrent;
  }

}
