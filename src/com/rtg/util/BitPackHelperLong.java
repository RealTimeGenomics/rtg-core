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
 * Helper for bit packing into long, only handles positive numbers
 */
public class BitPackHelperLong {

  /**
   * Number of bits in a long
   */
  public static final int LONG_BITS = Long.SIZE;

  private final int[] mFields;
  private final int[] mShifts;
  private final long[] mMasks;

  /**
   * Constructs the helper
   * @param sizeFields the size of each corresponding field
   * @throws IllegalArgumentException if total size &gt; 64
   */
  public BitPackHelperLong(int... sizeFields) {
    mFields = sizeFields.clone();
    int total = 0;
    for (int field : mFields) {
      total += field;
    }
    if (total > LONG_BITS) {
      throw new IllegalArgumentException("Cannot pack values into " + LONG_BITS + " bits when " + total + " bits are required");
    }
    mShifts = new int[mFields.length];
    //mShifts[0] = 0
    for (int i = 1; i < mShifts.length; ++i) {
      mShifts[i] = mShifts[i - 1] + mFields[i - 1];
    }
    mMasks = new long[mFields.length];
    for (int i = 0; i < mMasks.length; ++i) {
      mMasks[i] = (1L << mFields[i]) - 1;
    }
  }

  /**
   * Get the value of item in given field
   * @param fieldId id of field <code>0 &lt;= fieldId &lt; numberFields</code>
   * @param packedValue value of packed values
   * @return value of requested field
   */
  public long getField(int fieldId, long packedValue) {
    return (packedValue >> mShifts[fieldId]) & mMasks[fieldId];
  }

  /**
   * Pack all values into a long. NOTE: Do not try to store negative numbers
   * @param values values in order initially supplied to constructor
   * @return packed value
   */
  public long packValues(final long[] values) {
    assert values.length == mFields.length;
    long result = values[0];
    for (int i = 1; i < mShifts.length; ++i) {
      result |= values[i] << mShifts[i];
    }
    return result;
  }

}
