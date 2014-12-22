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


/**
 * Encapsulates parameters of a <code>NormalDistribution</code> over
 * the integers.
 *
 */
public class NormalDistribution {

  private long mCount;
  private double mSum;
  private double mSum2;

  /**
   * Create an empty normal distribution.
   */
  public NormalDistribution() { }

  /**
   * Create a normal distribution with provided values.
   * @param count number of observed values
   * @param sum sum of the observed values
   * @param sumSq sum of the observed values squared
   */
  public NormalDistribution(long count, double sum, double sumSq) {
    mCount = count;
    mSum = sum;
    mSum2 = sumSq;
  }

  /**
   * Merge a normal distribution into this one.
   * @param dist the other distribution
   */
  public void add(NormalDistribution dist) {
    mCount += dist.mCount;
    mSum += dist.mSum;
    mSum2 += dist.mSum2;
  }

  /**
   * Add a new data point to this distribution.
   * @param value the value
   */
  public void add(int value) {
    mCount++;
    mSum += value;
    mSum2 += value * value;
  }

  /**
   * @return the number of data values seen
   */
  public long count() {
    return mCount;
  }

  /**
   * @return the sum of the data values seen
   */
  public double sum() {
    return mSum;
  }

  /**
   * @return the sum of the data values squared
   */
  public double sumSq() {
    return mSum2;
  }

  /**
   * @return the mean of the normal distribution
   */
  public double mean() {
    if (mCount < 1) {
      throw new IllegalStateException("Count is not positive");
    }
    return mSum / mCount;
  }

  /**
   * @return the standard deviation of the normal distribution
   */
  public double stdDev() {
    if (mCount <= 1) {
      throw new IllegalStateException("Count is not greater than one");
    }
    return Math.sqrt((mSum2 - mSum * mean()) / (mCount - 1));
  }
}

