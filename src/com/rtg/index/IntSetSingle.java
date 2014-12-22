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
public final class IntSetSingle extends IntSet {

  private final int mRange;

  private final int mCapacity;

  private final int[] mValues;

  private final BitVector mBits;

  private int mNext = 0;

  /**
   * @param range of values that can be added.
   * @param capacity maximum number of values that can be added before <code>call</code>
   * is called on each add.
   * @param caller used to call action once or more for each unique value.
   */
  public IntSetSingle(final long range, final int capacity, final IntSetCaller caller) {
    super(caller);
    //System.err.println("IntSet range=" + range + " capacity=" + capacity);
    if (range > Integer.MAX_VALUE) {
      throw new RuntimeException("Range too large:" + range);
    }
    mRange = (int) range;
    mCapacity = capacity;
    mValues = new int[mCapacity];
    mBits = new BitVector(range);
  }

  @Override
  public void add(final int v) throws IOException {
    if (mBits.get(v)) {
      return;
    }
    if (mNext >= mCapacity) {
      call(v);
      return;
    }
    mValues[mNext] = v;
    mBits.set(v);
    mNext++;
  }

  @Override
  public void iterateClear() throws IOException {
    for (int i = 0; i < mNext; i++) {
      final int v = mValues[i];
      mBits.reset(v);
      call(v);
    }
    mNext = 0;
  }

  @Override
  public void iterateClearAll() throws IOException {
    iterateClear();
  }

  @Override
  public void toString(final StringBuilder sb) {
    sb.append("IntSet ").append(mNext).append(StringUtils.LS);
    for (int i = 0; i < mNext; i++) {
      sb.append(mValues[i]).append(" ");
    }
    sb.append(StringUtils.LS);
  }

  @Override
  public boolean globalIntegrity() {
    integrity();
    for (int i = 0; i < mNext; i++) {
      Exam.assertTrue(mBits.get(mValues[i]));
    }
    long total = 0;
    for (int i = 0; i < mRange; i++) {
      if (mBits.get(i)) {
        total++;
      }
    }
    Exam.assertTrue(total == mNext);
    return true;
  }

  @Override
  public boolean integrity() {
    Exam.assertTrue(mRange >= 0);
    Exam.assertTrue(mCapacity >= 0);
    Exam.assertTrue(mValues.length == mCapacity);
    Exam.assertTrue(mBits.length() == mRange);
    Exam.assertTrue(mNext >= 0 && mNext <= mCapacity);
    return true;
  }
}
