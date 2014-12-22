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
 * Bit pack two values into a long. NOTE: Do not try to store negative numbers
 */
public class BitPack2IntoLong {

  private final int mA;
  private final long mMaskA, mMaskB;

  /**
   * Constructs the helper
   * @param a size of first field.
   * @param b size of second field.
   * @throws IllegalArgumentException if total size &gt; 64
   */
  public BitPack2IntoLong(int a, int b) {
    if (a + b > 64) {
      throw new IllegalArgumentException();
    }
    mA = b;
    mMaskA = (1L << a) - 1;
    mMaskB = (1L << b) - 1;
  }

  /**
   * Get the value of item in given field.
   * @param fieldId id of field <code>0 &lt;= fieldId &lt; numberFields</code>
   * @param packedValue value of packed values
   * @return value of requested field
   */
  public long getField(final int fieldId, final long packedValue) {
    switch (fieldId) {
    case 0:
      return (packedValue >>> mA) & mMaskA;
    default:
      return packedValue & mMaskB;
    }
  }

  /**
   * Pack all values into a long. NOTE: Do not try to store negative numbers
   * @param a first value
   * @param b second value
   * @return packed value
   */
  public long packValues(final long a, final long b) {
    return (a << mA) | b;
  }

}
