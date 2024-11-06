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


import java.util.List;

/**
 * Encapsulates parameters of a <code>NormalDistribution</code> over the integers.
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
   * Combine a number of distributions into one
   * @param dists the distributions to combine
   */
  public void add(List<NormalDistribution> dists) {
    for (final NormalDistribution deviation : dists) {
      add(deviation);
    }
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
  public void add(long value) {
    ++mCount;
    mSum += value;
    mSum2 += Math.pow(value, 2);
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
      return 0.0;
    }
    return mSum / mCount;
  }

  /**
   * @return the standard deviation of the normal distribution
   */
  public double stdDev() {
    if (mCount <= 1) {
      return 0.0;
    }
    return Math.sqrt((mSum2 - mSum * mean()) / (mCount - 1));
  }

  @Override
  public String toString() {
    return mCount + "\t" + mean() + "\t" + stdDev() + "\t" + mSum + "\t" + mSum2;
  }

}

