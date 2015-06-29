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

import static com.rtg.util.StringUtils.LS;

/**
 * Class for calculating means, standard deviations and other simple statistics.
 */
public class Statistics {
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
      for (int i = 0; i < histogram.getLength(); i++) {
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
    mCount++;

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
