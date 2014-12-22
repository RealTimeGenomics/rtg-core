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

import java.io.Serializable;

/**
 * Class to hold either an absolute value or a percentage
 *
 */
public class IntegerOrPercentage implements Serializable, Comparable<IntegerOrPercentage> {
  private final int mValue;
  private final boolean mIsPercentage;

  /**
   * Constructs a new IntegerOrPercentage based on the supplied String.
   * If the string ends with "%", it is a percentage, otherwise an absolute value.
   *
   * @param s the string to parse
   * @throws NumberFormatException if not a valid number
   * @return the new instance
   */
  public static IntegerOrPercentage valueOf(final String s) {
    return new IntegerOrPercentage(s);
  }

  /**
   * Constructs a new IntegerOrPercentage based on the supplied int.
   *
   * @param i the int
   * @return the new instance
   */
  public static IntegerOrPercentage valueOf(final int i) {
    return new IntegerOrPercentage(i);
  }

  /**
   * Constructs a new IntegerOrPercentage based on the supplied String.
   * If the string ends with "%", it is a percentage, otherwise an absolute value.
   *
   * @param s the string to parse
   * @throws NumberFormatException if not a valid number
   */
  public IntegerOrPercentage(final String s) {
    final int pos = s.indexOf('%');
    if (pos >= 0) {
      mValue = Integer.parseInt(s.substring(0, pos));
      mIsPercentage = true;
    } else {
      mValue = Integer.parseInt(s);
      mIsPercentage = false;
    }
  }

  /**
   * Constructs a new IntegerOrPercentage based on the supplied int.
   *
   * @param i the int
   */
  public IntegerOrPercentage(final int i) {
    mValue = i;
    mIsPercentage = false;
  }

  /**
   * Computes the value of this IntegerOrPercentage.
   *
   * @param size the 100% value
   * @return the value unchanged if not a percentage, otherwise the percentage of size
   */
  public int getValue(final int size) {
    if (mIsPercentage) {
      return mValue * size / 100;
    } else {
      return mValue;
    }
  }

  /**
   * @return the raw value the user specified, independent of whether it is a percentage.
   */
  public int getRawValue() {
    return mValue;
  }

  /**
   *
   * @return if this is a percentage value
   */
  public boolean isPercentage() {
    return mIsPercentage;
  }

  /**
   * Compares with another IntegerOrPercentage
   * @param o an IntegerOrPercentage
   * @return &lt;0/0/&gt;0 if less than/equal/greater than
   */
  @Override
  public int compareTo(IntegerOrPercentage o) {
    if (mIsPercentage == o.mIsPercentage) {
      return mValue - o.mValue;
    } else {
      return mIsPercentage ? -1 : 1;
    }
  }

  @Override
  public boolean equals(Object o) {
    return o instanceof IntegerOrPercentage && compareTo((IntegerOrPercentage) o) == 0;
  }

  @Override
  public int hashCode() {
    int hash = 7;
    hash = 31 * hash + this.mValue;
    hash = 31 * hash + (this.mIsPercentage ? 1 : 0);
    return hash;
  }

  @Override
  public String toString() {
    return mIsPercentage ? Integer.toString(mValue) + "%" : Integer.toString(mValue);
  }
}
