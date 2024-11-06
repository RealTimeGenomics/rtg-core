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

import java.io.IOException;

import com.rtg.util.StringUtils;
import com.rtg.util.integrity.Exam;

/**
 * Allows storage and retrieval of a non-redundant set of integers.
 */
public final class IntSetWindow extends IntSet {

  private final int mRange;

  private final int mCapacity;

  private final int mWindow;

  private final int[][] mValues;

  private final BitVector mBits;

  private int mFirst;

  private int mLast;

  private final int[] mNext;

  /**
   * @param range of values that can be added.
   * @param capacity maximum number of values that can be added before <code>call</code>
   * @param window length of sliding window to use.
   * @param caller used to do calls for each unique value.
   * is called on each add.
   */
  public IntSetWindow(final long range, final int capacity, final int window, final IntSetCaller caller) {
    super(caller);
    //System.err.println("IntSet range=" + range + " capacity=" + capacity);
    if (range > Integer.MAX_VALUE) {
      throw new RuntimeException("Range too large:" + range);
    }
    mRange = (int) range;
    mCapacity = capacity;
    mWindow = window;
    mValues = new int[mWindow][mCapacity];
    mNext = new int[mWindow];
    mFirst = 0;
    mLast = mWindow - 1;
    mBits = new BitVector(range);
  }

  @Override
  public void add(final int v) throws IOException {
    if (mBits.get(v)) {
      return;
    }
    if (mNext[mLast] >= mCapacity) {
      call(v);
      return;
    }
    mValues[mLast][mNext[mLast]] = v;
    mBits.set(v);
    mNext[mLast]++;
  }

  @Override
  public void iterateClear() throws IOException {
    final int[] values = mValues[mFirst];
    for (int i = 0; i < mNext[mFirst]; ++i) {
      final int v = values[i];
      mBits.reset(v);
      call(v);
    }
    mNext[mFirst] = 0;
    ++mFirst;
    if (mFirst == mWindow) {
      mFirst = 0;
    }
    ++mLast;
    if (mLast == mWindow) {
      mLast = 0;
    }
  }

  @Override
  public void iterateClearAll() throws IOException {
    for (int i = 0; i < mWindow; ++i) {
      iterateClear();
    }
  }

  @Override
  public void toString(final StringBuilder sb) {
    sb.append("IntSetWindow ").append(mWindow).append(StringUtils.LS);
    for (int i = 0; i < mWindow; ++i) {
      final int j = mFirst + i;
      final int k = j >= mWindow ? j - mWindow : j;
      final int next = mNext[k];
      sb.append(next).append(": ");
      for (int l = 0; l < next; ++l) {
        sb.append(mValues[k][l]).append(" ");
      }
      sb.append(LS);
    }
  }

  @Override
  public boolean globalIntegrity() {
    integrity();
    for (int i = 0; i < mWindow; ++i) {
      for (int j = 0; j < mNext[i]; ++j) {
        Exam.assertTrue(mBits.get(mValues[i][j]));
      }
    }
    long totalv = 0;
    for (int i = 0; i < mRange; ++i) {
      if (mBits.get(i)) {
        ++totalv;
      }
    }
    long totalr = 0;
    for (int i = 0; i < mWindow; ++i) {
      totalr += mNext[i];
    }
    Exam.assertTrue(totalv == totalr);
    return true;
  }

  @Override
  public boolean integrity() {
    Exam.assertTrue(mRange >= 0);
    Exam.assertTrue(mCapacity >= 0);
    Exam.assertTrue(mValues.length == mWindow);
    for (int i = 0; i < mWindow; ++i) {
      Exam.assertTrue(mNext[i] >= 0 && mNext[i] <= mCapacity);
      Exam.assertTrue(mValues[i].length == mCapacity);
    }
    Exam.assertTrue(mBits.length() == mRange);
    return true;
  }
}
