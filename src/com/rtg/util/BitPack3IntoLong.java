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
 * Bit pack four values into a long. NOTE: Do not try to store negative numbers
 */
public class BitPack3IntoLong {

  private final int mA, mB ;
  private final long mMaskA, mMaskB, mMaskC;

  /**
   * Constructs the helper
   * @param a size of first field.
   * @param b size of second field.
   * @param c size of third field.
   * @throws IllegalArgumentException if total size &gt; 64
   */
  public BitPack3IntoLong(int a, int b, int c) {
    if (a + b + c > 64) {
      throw new IllegalArgumentException();
    }
    mB = c;
    mA = mB + b;
    mMaskA = (1L << a) - 1;
    mMaskB = (1L << b) - 1;
    mMaskC = (1L << c) - 1;
  }

  /**
   * Get the value of item in given field
   * @param fieldId id of field <code>0 &lt;= fieldId &lt; numberFields</code>
   * @param packedValue value of packed values
   * @return value of requested field
   */
  public long getField(final int fieldId, final long packedValue) {
    switch (fieldId) {
    case 0:
      return (packedValue >>> mA) & mMaskA;
    case 1:
      return (packedValue >>> mB) & mMaskB;
    default:
      return packedValue & mMaskC;
    }
  }

  /**
   * Pack all values into a long. NOTE: Do not try to store negative numbers
   * @param a first value
   * @param b second value
   * @param c third value
   * @return packed value
   */
  public long packValues(final long a, final long b, final long c) {
    return (a << mA) | (b << mB) | c;
  }

}
