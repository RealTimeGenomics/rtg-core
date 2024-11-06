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
 * Bit pack two values into a long. NOTE: Do not try to store negative numbers.
 */
public class BitPack2IntoLong {

  private final int mBitsForB;
  private final long mMaskB;

  /**
   * Constructs the helper
   * @param a size of first field.
   * @param b size of second field.
   * @throws IllegalArgumentException if total size &gt; 64
   */
  public BitPack2IntoLong(int a, int b) {
    if (a + b > 64 || a < 1 || b < 1) {
      throw new IllegalArgumentException();
    }
    mBitsForB = b;
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
      return packedValue >>> mBitsForB;
    default:
      return packedValue & mMaskB;
    }
  }

  /**
   * Pack all values into a long. NOTE: Do not try to store negative numbers.
   * @param a first value
   * @param b second value
   * @return packed value
   */
  public long packValues(final long a, final long b) {
    assert a >= 0 && b >= 0;
    return (a << mBitsForB) | b;
  }

}
