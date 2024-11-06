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
package com.rtg.index;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.util.StringUtils;

/**
 * A vector of bits.
 */
@TestClass(value = {"com.rtg.index.BitVectorTest"})
public abstract class AbstractBitVector {

  static final int BITS_PER_ENTRY = 5;

  static final int BITS_PER_BYTE = 3;

  static final int BYTES_PER_ENTRY = 1 << (BITS_PER_ENTRY - BITS_PER_BYTE);

  static final int MASK = (1 << BITS_PER_ENTRY) - 1;

  /**
   * Create a new bit vector choosing the best implementation depending on the length.
   * @param length in bits of the vector.
   * @return the <code>AbstractBitVector</code>
   */
  public static AbstractBitVector createBitVector(final long length) {
    final long allocLength = (length + MASK) >> BITS_PER_ENTRY;
    if (allocLength >= (Integer.MAX_VALUE >> 2)) {
      return new BigBitVector(length);
    } else {
      return new BitVector(length);
    }
  }

  protected final long mLength;

  /**
   * Creates a new <code>BitVector</code> with the specified length.
   *
   * @param length the number of bits.
   */
  public AbstractBitVector(final long length) {
    if (length < 0) {
      throw new RuntimeException("Negative length:" + length);
    }
    mLength = length;
  }


  /**
   * Returns true if the specified bit is set.
   *
   * @param index the index of the bit.
   * @return true if the specified bit is set.
   */
  public abstract boolean get(long index);


  /**
   * Sets a bit.
   *
   * @param index the index of the bit.
   */
  public abstract void set(final long index);

  /**
   * Resets a bit (to 0).
   *
   * @param index the index of the bit.
   */
  public abstract void reset(long index);


  /**
   * Gets the number of bytes in the vector.
   *
   * @return number of bytes in bit vector - just takes account of the array.
   */
  public abstract long bytes();


  /**
   * Gets the length of the bit vector.
   *
   * @return number of bits.
   */
  public long length() {
    return mLength;
  }

  @Override
  public String toString() {
    final StringBuilder sb = new StringBuilder();
    sb.append("BitVector[").append(mLength).append("]").append(StringUtils.LS);
    for (long l = 0; l < mLength; l += 100) {
      sb.append("[").append(l).append("]\t");
      final long min = Math.min(l + 100, mLength);
      for (long m = l; m < min; ++m) {
        if (m % 10 == 0) {
          sb.append(" ");
        }
        sb.append(get(m) ? "1" : "0");
      }
      sb.append(StringUtils.LS);
    }
    sb.append(StringUtils.LS);
    return sb.toString();
  }
}
