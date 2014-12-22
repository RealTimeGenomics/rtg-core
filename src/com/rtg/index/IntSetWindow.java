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
    mValues = new int[mWindow][];
    for (int i = 0; i < mWindow; i++) {
      mValues[i] = new int[mCapacity];
    }
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
    for (int i = 0; i < mNext[mFirst]; i++) {
      final int v = values[i];
      mBits.reset(v);
      call(v);
    }
    mNext[mFirst] = 0;
    mFirst++;
    if (mFirst == mWindow) {
      mFirst = 0;
    }
    mLast++;
    if (mLast == mWindow) {
      mLast = 0;
    }
  }

  @Override
  public void iterateClearAll() throws IOException {
    for (int i = 0; i < mWindow; i++) {
      iterateClear();
    }
  }

  @Override
  public void toString(final StringBuilder sb) {
    sb.append("IntSetWindow ").append(mWindow).append(StringUtils.LS);
    for (int i = 0; i < mWindow; i++) {
      final int j = mFirst + i;
      final int k = j >= mWindow ? j - mWindow : j;
      final int next = mNext[k];
      sb.append(next).append(": ");
      for (int l = 0; l < next; l++) {
        sb.append(mValues[k][l]).append(" ");
      }
      sb.append(LS);
    }
  }

  @Override
  public boolean globalIntegrity() {
    integrity();
    for (int i = 0; i < mWindow; i++) {
      for (int j = 0; j < mNext[i]; j++) {
        Exam.assertTrue(mBits.get(mValues[i][j]));
      }
    }
    long totalv = 0;
    for (int i = 0; i < mRange; i++) {
      if (mBits.get(i)) {
        totalv++;
      }
    }
    long totalr = 0;
    for (int i = 0; i < mWindow; i++) {
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
    for (int i = 0; i < mWindow; i++) {
      Exam.assertTrue(mNext[i] >= 0 && mNext[i] <= mCapacity);
      Exam.assertTrue(mValues[i].length == mCapacity);
    }
    Exam.assertTrue(mBits.length() == mRange);
    return true;
  }
}
