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
