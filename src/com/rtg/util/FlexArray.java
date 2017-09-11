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

import java.util.Arrays;

/**
 * Expands array to required size on set.
 * @param <T> type of objects stored in array.
 */
public class FlexArray<T> {
  private Object[] mList;
  private int mSize;

  /**
   * Constructor.
   */
  public FlexArray() {
    this(10);
  }

  /**
   * @param initialCapacity initial capacity of internal array
   */
  public FlexArray(int initialCapacity) {
    assert initialCapacity >= 0;
    mList = new Object[initialCapacity];
  }

  static int newSize(int currentLength, int cap) {
    int newSize = currentLength + 2;
    while (newSize < cap) {
      final long newSizeL = ((long) newSize * 3L) / 2L;
      if (newSizeL > Integer.MAX_VALUE) {
        newSize = Integer.MAX_VALUE;
      } else {
        newSize = (int) newSizeL;
      }
    }
    return newSize;
  }

  /**
   * Set the arrays value at the given position, expanding the size and capacity if required.
   * @param index position to set
   * @param obj the value
   */
  public void set(int index, T obj) {
    ensureCapacity(index + 1);
    mList[index] = obj;
    if (mSize <= index) {
      mSize = index + 1;
    }
  }

  /**
   * @return the maximum index set + 1
   */
  public int size() {
    return mSize;
  }

  /**
   * Get an array of length <code>size()</code> containing the values
   * @param a array of same type to place the value in, if not the right length then a new one will be used instead
   * @return the array the values were copied to.
   */
  @SuppressWarnings("unchecked")
  public T[] toArray(T[] a) {
    if (a.length != mSize) {
      return (T[]) Arrays.copyOf(mList, mSize, a.getClass());
    }
    System.arraycopy(mList, 0, a, 0, mSize);
    return a;
  }

  private void ensureCapacity(int cap) {
    if (cap <= mList.length) {
      return;
    }
    final int newSize = newSize(mList.length, cap);
    mList = Arrays.copyOf(mList, newSize);
  }

  int capacity() {
    return mList.length;
  }
}
