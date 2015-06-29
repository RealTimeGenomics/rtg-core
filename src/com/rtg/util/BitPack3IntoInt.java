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
 * Bit pack three values into an integer. NOTE: Do not try to store negative numbers.
 */
public class BitPack3IntoInt {

  private final int mA, mB;
  private final int mMaskB, mMaskC;

  /**
   * Constructs the helper.
   * @param a size of first field.
   * @param b size of second field.
   * @param c size of third field.
   * @throws IllegalArgumentException if total size &gt; 32
   */
  public BitPack3IntoInt(final int a, final int b, final int c) {
    if (a + b + c > 32 || a < 1 || b < 1 || c < 1) {
      throw new IllegalArgumentException();
    }
    mB = c;
    mA = mB + b;
    mMaskB = (1 << b) - 1;
    mMaskC = (1 << c) - 1;
  }

  /**
   * Get the value of item in given field
   * @param fieldId id of field <code>0 &lt;= fieldId &lt; numberFields</code>
   * @param packedValue value of packed values
   * @return value of requested field
   */
  public int getField(final int fieldId, final int packedValue) {
    switch (fieldId) {
    case 0:
      return packedValue >>> mA;
    case 1:
      return (packedValue >>> mB) & mMaskB;
    default:
      return packedValue & mMaskC;
    }
  }

  /**
   * Pack all values into an integer.
   * It is the callers responsibility to ensure that the values do not exceed allowed
   * number of bits. NOTE: Do not try to store negative numbers.
   * @param a first value
   * @param b second value
   * @param c third value
   * @return packed value
   */
  public int packValues(final int a, final int b, final int c) {
    assert a >= 0 && b >= 0 && c >= 0;
    return (a << mA) | (b << mB) | c;
  }

}
