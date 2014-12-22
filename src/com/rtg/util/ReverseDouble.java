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
 * A Double where the natural ordering is reversed (e.g. for sorting).
 */
public class ReverseDouble implements Comparable<ReverseDouble> {

  private final double mValue;

  /**
   * @param value value to reverse
   */
  public ReverseDouble(final double value) {
    mValue = value;
    if (Double.isNaN(value)) {
      throw new IllegalArgumentException();
    }
  }

  @Override
  public int compareTo(final ReverseDouble that) {
    return Double.compare(that.mValue, this.mValue);
  }

  /**
   * Returns the double value.
   * @return the double value.
   */
  public double doubleValue() {
    return mValue;
  }

  @Override
  public boolean equals(final Object obj) {
    return (obj instanceof ReverseDouble)
        && Double.doubleToLongBits(((ReverseDouble) obj).mValue) == Double.doubleToLongBits(this.mValue);
  }

  @Override
  public int hashCode() {
    final long lb = Double.doubleToRawLongBits(mValue);
    return ((int) lb) ^ ((int) (lb >> 32));
  }

  @Override
  public String toString() {
    return Double.toString(mValue);
  }
}
