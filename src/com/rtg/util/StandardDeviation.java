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

import java.util.List;

/**
 * Calculates mean and standard deviation from samples.
 */
public class StandardDeviation {
  private double mSum = 0;
  private double mSumSquared = 0;
  private long mNumSamples = 0;

  /**
   * @param value the value of the new sample
   */
  public void addSample(long value) {
    mNumSamples++;
    mSum += value;
    mSumSquared += Math.pow(value, 2);
  }

  /**
   * Calculates standard deviation
   * <code> &#0963;<sup>2</sup> = (&#0931;(value<sub>i</sub><sup>2</sup>) - (&#0931;(value<sub>i</sub>)<sup>2</sup> / i)) / i</code>
   * @return standard deviation
   */
  public double standardDeviation() {
    if (mNumSamples == 0) {
      return 0.0;
    }
    final double term2 = Math.pow(mSum, 2) / mNumSamples;
    final double sqrt = Math.sqrt((mSumSquared - term2) / mNumSamples);
    if (Double.isNaN(sqrt)) {
      System.err.println("mSum=" + mSum
          + " mNumSamples=" + mNumSamples
          + " mSumSquared=" + mSumSquared
          + " term2=" + term2
          + " (mSumSquared - term2)=" + (mSumSquared - term2)
          + " ((mSumSquared - term2) / mNumSamples)=" + ((mSumSquared - term2) / mNumSamples)
      );
    }
    return sqrt;
  }

  /**
   * Calculates mean
   * @return mean
   */
  public double mean() {
    if (mNumSamples == 0) {
      return 0.0;
    }
    return mSum / mNumSamples;
  }

  @Override
  public String toString() {
    return mNumSamples + "\t" + mean() + "\t" + standardDeviation() + "\t" + mSum + "\t" + mSumSquared;
  }

  /**
   * Combine a number of distributions into one
   * @param deviations the distributions to combine
   * @return the sum of the distributions
   */
  public static StandardDeviation combine(List<StandardDeviation> deviations) {
    final StandardDeviation combined = new StandardDeviation();
    for (StandardDeviation deviation : deviations) {
      combined.mSum += deviation.mSum;
      combined.mSumSquared += deviation.mSumSquared;
      combined.mNumSamples += deviation.mNumSamples;
    }
    return combined;
  }

}
