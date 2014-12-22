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
 * Contains a value and a lower+higher bound pair
 */
public class BoundedDouble {
  double mValue;
  double mLow;
  double mHigh;

  /**
   * Construct a double with a low and high bounds
   * @param val the value
   * @param low the low bound
   * @param high the high bound
   */
  public BoundedDouble(double val, double low, double high) {
    assert low <= val && val <= high : low + "<=" + val + "<=" + high;
    mValue = val;
    mLow = low;
    mHigh = high;
  }

  public double getValue() {
    return mValue;
  }

  public double getLow() {
    return mLow;
  }

  public double getHigh() {
    return mHigh;
  }

  @Override
  public String toString() {
    return String.valueOf(mValue);
  }
}
