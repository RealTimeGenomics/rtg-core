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

package com.rtg.variant.sv;

import java.util.Arrays;

import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 * An SamCounts that supports efficient calculation of the sum of counts between two points.
 *
 */
public class CumulativeSamCounts extends IntegralAbstract implements SamCounts {

  private double[] mCounts;
  private int mLength;
  private int mOffset = 0; // Minimum valid external coordinate that can be accessed.
  private int mLastIncrement = 1; // Lowest internal coordinate that can be incremented (i.e. external coord + 1)

  private double[] mSumLn;
  private int mLnLastUpdated = 0;

  private final Corrections mCorrections;

  /**
   * Construct a SamArray
   * @param length the size of the array.
   * @param corrections optional coverage corrections object
   */
  public CumulativeSamCounts(int length, Corrections corrections) {
    mCounts = new double[length + 1];
    mSumLn = new double[length + 1];
    mLength = length;
    mCorrections = corrections;
  }

  @Override
  public void increment(final int index) {
    if (index < 0 || index >= length()) {
      throw new IndexOutOfBoundsException("index=" + index + " length=" + length());
    }
    if (index < mOffset || index >= (mOffset + mCounts.length)) {
      throw new IndexOutOfBoundsException("index=" + index + " offset=" + mOffset + " length=" + mCounts.length);
    }
    if (index < mLastIncrement - 1) {
      throw new IllegalStateException("Out of order incrementing not supported");
    }
    copyTo(index);
    final double incr;
    if (mCorrections == null) {
      incr = 1.0;
    } else {
      final double correction = mCorrections.correction(index);
      if (correction == 0.0) {
        incr = 0.0;
      } else {
        incr = 1.0 / correction;
      }
    }
    mCounts[mLastIncrement % mCounts.length] += incr;
  }

  private void copyTo(int index) {
    final double tot = mCounts[mLastIncrement % mCounts.length];
    while (mLastIncrement <= index)   {
      ++mLastIncrement;
      mCounts[mLastIncrement % mCounts.length] = tot;
    }
  }

  @Override
  public double count(int base, int index) {
    if (base + index < mOffset) {
      return 0;
    }
    return count(base, index, index + 1);
  }

  @Override
  public double count(int base, int index, int index2) {
    assert index2 > index;
    final int ex1 = Math.min(Math.max(base + index, mOffset), mLastIncrement);
    final int ex2 = Math.min(Math.max(base + index2, mOffset), mLastIncrement);
    if ((ex1 < mOffset) || (ex2 < mOffset)) {
      throw new IllegalArgumentException("Accessing coordinate out of range min=" + mOffset + " base=" + base + " index=" + index + " index2=" + index2);
      //return 0;
    }
    // count at ex1 = counts[ex1 + 1] - counts[ex1];
    // count at ex2 = counts[ex2 + 1] - counts[ex2];
    // cumulative to ex1 noninclusive = counts[ex1];
    // cumulative to ex2 noninclusive = counts[ex2];
    // cumulative from ex1 inclusive to ex2 noninclusive = counts[ex2] - counts[ex1]
    return mCounts[ex2 % mCounts.length] - mCounts[ex1 % mCounts.length];
  }


  @Override
  public double sumLn(int base, int index, int index2) {
    assert index2 > index;

    // Update calculation from last point
    while (mLnLastUpdated < (mLastIncrement - 1)) {
      ++mLnLastUpdated;
      updateLogTerm(mLnLastUpdated);
    }
    updateLogTerm(mLastIncrement); // Recalculate at this point in case more has been added since last time

    final int ex1 = Math.min(Math.max(base + index, mOffset), mLastIncrement);
    final int ex2 = Math.min(Math.max(base + index2, mOffset), mLastIncrement);
    if ((ex1 < mOffset) || (ex2 < mOffset)) {
      throw new IllegalArgumentException("Accessing coordinate out of range min=" + mOffset + " base=" + base + " index=" + index + " index2=" + index2);
    }
    return mSumLn[ex2 % mSumLn.length] - mSumLn[ex1 % mSumLn.length];
  }

  private void updateLogTerm(int index) {
    final int i = index % mCounts.length;
    final int im1 = (index - 1) % mCounts.length;
    final double ni = mCounts[i] - mCounts[im1];
    mSumLn[i] = mSumLn[im1] + ((ni > 0) ? ni * Math.log(ni) : 0);
  }

  @Override
  public void reset(final int length, int bufferSize) {
    if (mCounts.length < bufferSize + 1) {
      mCounts = new double[bufferSize + 1];
      mSumLn = new double[bufferSize + 1];
    } else {
      Arrays.fill(mCounts, 0);
      Arrays.fill(mSumLn, 0);
    }
    mLength = length;
    mOffset = 0;
    mLastIncrement = 1;
    mLnLastUpdated = 0;
  }

  @Override
  public void flushTo(int offset) {
    if (offset < mOffset) {
      throw new IllegalArgumentException("Attempt to flush to position earlier than already flushed: " + offset + " < " + mOffset);
    }
    if (offset - mOffset >= mCounts.length) {
      Arrays.fill(mCounts, 0);
      Arrays.fill(mSumLn, 0);
    } else {
      final int ix1 = mOffset % mCounts.length;
      final int ix2 = offset % mCounts.length;
      if (ix2 >= ix1) {
        Arrays.fill(mCounts, ix1, ix2, 0);
        Arrays.fill(mSumLn, ix1, ix2, 0);
      } else {
        Arrays.fill(mCounts, ix1, mCounts.length, 0);
        Arrays.fill(mCounts, 0, ix2, 0);
        Arrays.fill(mSumLn, ix1, mSumLn.length, 0);
        Arrays.fill(mSumLn, 0, ix2, 0);
      }
    }
    mOffset = offset;
    mLastIncrement = Math.max(mOffset, mLastIncrement);
    mLnLastUpdated = Math.max(mOffset, mLnLastUpdated);
  }

  @Override
  public int length() {
    return mLength;
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
