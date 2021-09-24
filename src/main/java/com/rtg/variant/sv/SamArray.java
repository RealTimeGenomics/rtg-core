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

package com.rtg.variant.sv;

import java.util.Arrays;

import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 * An array of counts.
 *
 */
public class SamArray extends IntegralAbstract implements SamCounts {

  private float[] mCounts;
  private int mLength;
  private int mOffset = 0;

  /**
   * Construct a SamArray
   * @param length the size of the array.
   */
  public SamArray(int length) {
    mCounts = new float[length];
    mLength = length;
  }

  @Override
  public void increment(final int index) {
    increment(index, 1.0);
  }

  /**
   * Increment the count at a position in the array.
   * @param index position used.
   * @param incr the amount to increment by
   */
  protected void increment(final int index, final double incr) {
    if (index < 0 || index >= length()) {
      throw new IndexOutOfBoundsException("index=" + index + " length=" + length());
    }
    if (index < mOffset || index >= (mOffset + mCounts.length)) {
      throw new IndexOutOfBoundsException("index=" + index + " offset=" + mOffset + " length=" + mCounts.length);
    }
    final int index2 = index % mCounts.length;
    mCounts[index2] += (float) incr;
  }

  @Override
  public double count(int base, int index) {
    final int ix = base + index;
    if (ix < mOffset || ix >= (mOffset + mCounts.length)) {
      return 0.0;
    }
    final int ix2 = ix % mCounts.length;
    return mCounts[ix2];
  }

  @Override
  public double count(int base, int index, int index2) {
    double sum = 0;
    for (int i = index; i < index2; ++i) {
      sum += count(base, i);
    }
    return sum;
  }

  @Override
  public double sumLn(int base, int index, int index2) {
    double sum = 0;
    for (int i = index; i < index2; ++i) {
      final double count = count(base, i);
      if (count > 0) {
        sum += count * Math.log(count);
      }
    }
    return sum;
  }

  @Override
  public void reset(final int length, int bufferSize) {
    if (mCounts.length < bufferSize) {
      mCounts = new float[bufferSize];
    } else {
      Arrays.fill(mCounts, 0);
    }
    mLength = length;
    mOffset = 0;
  }

  @Override
  public void flushTo(int offset) {
    if (offset < mOffset) {
      throw new IllegalArgumentException("Attempt to flush to position earlier than already flushed");
    }
    if (offset - mOffset >= mCounts.length) {
      Arrays.fill(mCounts, 0);
    } else {
      final int ix1 = mOffset % mCounts.length;
      final int ix2 = offset % mCounts.length;
      if (ix2 >= ix1) {
        Arrays.fill(mCounts, ix1, ix2, 0);
      } else {
        Arrays.fill(mCounts, ix1, mCounts.length, 0);
        Arrays.fill(mCounts, 0, ix2, 0);
      }
    }
    mOffset = offset;
  }

  @Override
  public int length() {
    return mLength;
  }

  /**
   * Produce a new array swapped around pivot and then offset.
   * The values at the start or end are extended to ensure that
   * the whole of the new array has been filled.
   * Useful for testing.
   * @param pivot the point around which swapping is done.
   * @param offset to be subtracted after swapping.
   * @return the reversed array.
   */
  SamArray reverse(final int pivot, final int offset) {
    final SamArray rev = new SamArray(mLength);
    final int c = 2 * pivot - offset - 1;
    //extend at top of new array if necessary
    final double lov = count(0, 0);
    for (int i = c + 1; i < mLength; ++i) {
      rev.increment(i, lov);
    }
    //extend at bottom of new array if necessary
    final double hiv = count(0, mLength - 1);
    final int cc = c - mLength + 1;
    for (int i = 0; i < cc; ++i) {
      rev.increment(i, hiv);
    }
    //fill in the bulk of the array
    for (int i = 0; i < mLength; ++i) {
      final int j = c - i;
      if (j >= 0 && j < mLength) {
        rev.increment(j, count(0, i));
      }
    }
    return rev;
  }

  @Override
  public void toString(StringBuilder sb) {
    sb.append(Arrays.toString(mCounts));
  }

  @Override
  public boolean integrity() {
    Exam.assertTrue(mCounts.length >= length());
    return true;
  }
}
