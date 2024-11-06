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

import static com.rtg.util.StringUtils.LS;

/**
 * Class for calculating means, standard deviations and other simple statistics.
 */
public final class Statistics {
  private long mSum;
  private long mSumSquares;
  private long mCount;
  private long mMax;
  private long mMin;

  private Histogram mHistogram = null;

  /**
   * Create a statistics object resetting internal counters.
   */
  public Statistics() {
    this(null);
  }

  /**
   * Create statistics object with a backing histogram to collect more statistics.
   *
   * @param histogram backing structure to retain data
   */
  public Statistics(Histogram histogram) {
    reset();
    if (histogram != null) {
      mHistogram = histogram;
      for (int i = 0; i < histogram.getLength(); ++i) {
        final long count = histogram.getValue(i);
        if (count > 0) {
          addSample(i, count);
        }
      }
    }
  }

  /**
   * Reset statistics object, removing existing data.
   */
  public void reset() {
    mSum = 0L;
    mSumSquares = 0L;
    mCount = 0L;
    mMax = Long.MIN_VALUE;
    mMin = Long.MAX_VALUE;
    if (mHistogram != null) {
      mHistogram = new Histogram(); // only way to clear existing histogram...
    }
  }

  /**
   * Add a new value to the collection.
   *
   * @param x value to add
   */
  public void addSample(int x) {
    mSum += x;
    mSumSquares += (long) x * (long) x;
    mMin = x < mMin ? x : mMin;
    mMax = x > mMax ? x : mMax;
    ++mCount;

    // If the sum of squares exceeds Long.MAX_VALUE, this means the
    // value has overflowed; reset the state back to zero and start again.
    // All previous calculations are lost.  (Better as all doubles?)
    if (mSumSquares < 0L) {
      reset();
    }
  }

  /**
   * Add a new value to the collection a multiple number of times.
   *
   * @param x value to add
   * @param times number of times to add
   */
  public void addSample(int x, long times) {
    mSum += (long) x * times;
    mSumSquares += (long) x * (long) x * times;
    mMin = x < mMin ? x : mMin;
    mMax = x > mMax ? x : mMax;
    mCount += times;

    // If the sum of squares exceeds Long.MAX_VALUE, this means the
    // value has overflowed; reset the state back to zero and start again.
    // All previous calculations are lost.  (Better as all doubles?)
    if (mSumSquares < 0L) {
      reset();
    }
  }

  /**
   * Returns the mean of the numbers collected.
   *
   * @return mean value of collection
   */
  public double mean() {
    return mCount > 0L ? (double) mSum / mCount : 0.0;
  }

  /**
   * Returns the standard deviation of the numbers in the collection.
   *
   * @return standard deviation of collection
   */
  public double standardDeviation() {
    return Math.sqrt(variance());
  }

  /**
   * Returns the variance of the numbers in the collection.
   *
   * @return variance of collection
   */
  public double variance() {
    return mCount > 1L ? (mSumSquares - (double) mSum * mSum / mCount) / (mCount - 1) : 0.0;
  }


  /**
   * Returns the number of values in the collection.
   *
   * @return number of values
   */
  public long count() {
    return mCount;
  }

  /**
   * Returns the sum of the numbers in the collection.
   *
   * @return sum of collection
   */
  public long sum() {
    return mSum;
  }

  /**
   * Returns the maximum value in the collection.
   *
   * @return maximum value in collection
   */
  public long max() {
    return mCount > 0L ? mMax : 0L;
  }

  /**
   * Returns the minimum value in the collection.
   *
   * @return minimum value in collection
   */
  public long min() {
    return mCount > 0L ? mMin : 0L;
  }

  @Override
  public String toString()    {
    // var
    // se
    // skew
    // kurt

    // median
    // mode
    return "n    : " + mCount + LS + "min  : " + mMin + LS + "max  : " + mMax + LS + "sum  : " + mSum + LS + "ss   : " + mSumSquares + LS + "mean : " + Utils.realFormat(mean(), 4) + LS + "sd   : " + Utils.realFormat(standardDeviation(), 4) + LS;
  }
}
