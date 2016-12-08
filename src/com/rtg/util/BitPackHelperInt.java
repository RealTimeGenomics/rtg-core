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
 * Helper for bit packing into int only works with positive numbers
 */
public class BitPackHelperInt {

  /**
   * Number of bits in a int
   */
  public static final int INT_BITS = Integer.SIZE;

  private final int[] mFields;
  private final int[] mShifts;
  private final int[] mMasks;

  /**
   * Constructs the helper
   * @param sizeFields the size of each corresponding field
   * @throws IllegalArgumentException if total size &gt; 64
   */
  public BitPackHelperInt(int... sizeFields) {
    mFields = sizeFields.clone();
    int total = 0;
    for (int field : mFields) {
      total += field;
    }
    if (total > INT_BITS) {
      throw new IllegalArgumentException("Cannot pack values into " + INT_BITS + " bits when " + total + " bits are required");
    }
    mShifts = new int[mFields.length];
    //mShifts[0] = 0;
    for (int i = 1; i < mShifts.length; ++i) {
      mShifts[i] = mShifts[i - 1] + mFields[i - 1];
    }
    mMasks = new int[mFields.length];
    for (int i = 0; i < mMasks.length; ++i) {
      mMasks[i] = (1 << mFields[i]) - 1;
    }
  }

  /**
   * Get the value of item in given field
   * @param fieldId id of field <code>0 &lt;= fieldId &lt; numberFields</code>
   * @param packedValue value of packed values
   * @return value of requested field
   */
  public int getField(int fieldId, int packedValue) {
    return (packedValue >> mShifts[fieldId]) & mMasks[fieldId];
  }

  /**
   * Pack all values into a int. NOTE: Do not try to store negative numbers
   * @param values values in order initially supplied to constructor
   * @return packed value
   */
  public int packValues(final int[] values) {
    assert values.length == mFields.length;
    int result = values[0];
    for (int i = 1; i < mShifts.length; ++i) {
      result |= values[i] << mShifts[i];
    }
    return result;
  }

}
