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

package com.rtg.util.integrity;

/**
 * Used to check if integer values are in a (contiguous) integer subrange.
 * Includes methods for converting to and from strings and also allows
 * an "invalid" value. Expected usage is that a singleton will be created for each
 * range that will then be used in assert statements, as well as input and output conversion.
 * The invalid value is assumed (this is checked by the integrity method) to be one less than the smallest valid value.
 */
public class IntegerRange extends IntegralAbstract {

  private final int mLo;

  private final int mHi;

  private final boolean mHasInvalid;

  private final int mInvalid;

  protected IntegerRange(final int lo, final int hi) {
    this(false, Integer.MIN_VALUE, lo, hi);
  }

  protected IntegerRange(final int invalid, final int lo, final int hi) {
    this(true, invalid, lo, hi);
    integrity();
  }

  protected IntegerRange(final boolean hasInvalid, final int invalid, final int lo, final int hi) {
    mLo = lo;
    mHi = hi;
    mHasInvalid = hasInvalid;
    if (hasInvalid) {
      mInvalid = invalid;
    } else {
      mInvalid = Integer.MIN_VALUE;
    }
    integrity();
  }

  /**
   * Get the lowest value for any valid integer.
   * Does not include the invalid value in the range from <code>low()</code> to <code>high()</code>.
   * @return the lower bound (inclusive).
   */
  public int low() {
    return mLo;
  }

  /**
   * Get the highest value (inclusive) for any valid integer.
   * @return the upper bound (inclusive).
   * Does not include the invalid value in the range from <code>low()</code> to <code>high()</code>.
   */
  public int high() {
    return mHi;
  }

  /**
   * @return true iff an explicit invalid value is defined.
   */
  public boolean hasInvalid() {
    return mHasInvalid;
  }

  /**
   * @return the invalid value.
   * @exception RuntimeException if there is no invalid value defined.
   */
  public int invalid() {
    if (!mHasInvalid) {
      throw new RuntimeException("No invalid value defined");
    }
    return mInvalid;
  }

  /**
   * Check that i is in range (excluding the invalid value if any).
   * @param i value to be checked.
   * @return true iff if in range and not equal to invalid value.
   */
  public final boolean valid(final long i) {
    return mLo <= i && i <= mHi;
  }

  /**
   * Assert that i is in range and not equal to the invalid value (if any).
   * @param i value to be checked.
   * @return true so can be used in assert statements.
   */
  public final boolean checkValid(final int i) {
    Exam.assertTrue(i + ":" + toString(), valid(i));
    return true;
  }

  /**
   * Check if i is in range (includes invalid value).
   * @param i value to be checked.
   * @return true iff i is in range.
   */
  public final boolean inRange(final long i) {
    return valid(i) || (hasInvalid() && i == mInvalid);
  }

  /**
   * Assert that i is in range including invalid value.
   * @param i value to be checked.
   * @return true so can be used in assert statements.
   */
  public final boolean check(final long i) {
    Exam.assertTrue(i + ":" + toString(), inRange(i));
    return true;
  }

  /**
   * Convert i to a string.
   * Individual implementations may well want to override this.
   * @param i value to be converted.
   * @return string representation of i.
   */
  public String toString(final int i) {
    return "" + i;
  }

  /**
   * Get value of str as an integer.
   * Should obey the equation <code>i = valueOf(toString(i))</code>.
   * @param str to be converted.
   * @return value.
   */
  public int valueOf(final String str) {
    final int i = Integer.parseInt(str);
    check(i);
    return i;
  }

  @Override
  public void toString(final StringBuilder sb) {

    if (mHasInvalid) {
      sb.append(mInvalid).append("i");
    }
    sb.append("[").append(mLo).append("...").append(mHi).append("]");
  }

  @Override
  public boolean integrity() {
    Exam.assertTrue(mLo <= mHi);
    if (mHasInvalid) {
      Exam.assertTrue(mLo - 1 == mInvalid);
    } else {
      Exam.assertEquals(Integer.MIN_VALUE, mInvalid);
    }
    return true;
  }

}
