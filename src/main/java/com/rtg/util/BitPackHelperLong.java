/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
