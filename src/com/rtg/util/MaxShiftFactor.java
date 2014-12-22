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
 * Works out maximum aligner band width based on read length and a scaling factor.
 *
 */
public class MaxShiftFactor {
  private final double mFactor;

  /**
   * Set scaling factor.
   * @param factor scaling factor.
   */
  public MaxShiftFactor(double factor) {
    if (factor < 0.0 || factor > 1) {
      throw new IllegalArgumentException("factor must be between 0 and 1: " + factor);
    }
    mFactor = factor;
  }

  /**
   * Returns the maximum shift for the give alignment threshold.
   * @param readLength the read length, to use for max shift threshold
   * @return maximum shift that is needed
   */
  public int calculateMaxShift(int readLength) {
    return Math.max((int) (readLength * mFactor), 7);
  }

  /**
   * Return the shift factor.
   * @return max shift factor.
   */
  public double getFactor() {
    return mFactor;
  }

  /**
   */
  @Override
  public boolean equals(Object o) {
    return o != null && o instanceof MaxShiftFactor && Double.doubleToLongBits(((MaxShiftFactor) o).mFactor) == Double.doubleToLongBits(mFactor);
  }

  /**
   */
  @Override
  public int hashCode() {
    return Double.valueOf(mFactor).hashCode() - 7;
  }

}
