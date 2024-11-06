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

import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 */
public final class HashBitHandle extends IntegralAbstract {

  private final int mHashBits;

  private final int mLengthBits;


  /**
   * @param hashBits the number of valid bits used in the original hash.
   * @param lengthBits the number of bits to be used to address the vector.
   */
  public HashBitHandle(final int hashBits, final int lengthBits) {
    mHashBits = hashBits;
    mLengthBits = lengthBits;
    integrity();
  }

  /**
   * Create a <code>HashBitVector</code> with the specified hash and length.
   * @return a <code>HashBitVector</code> with the specified hash and length.
   */
  public HashBitVector create() {
    return new HashBitVector(mHashBits, mLengthBits);
  }

  /**
   * Get the length of the <code>HashBitVector</code> in bits.
   * @return the length of the <code>HashBitVector</code> in bits.
   */
  public long length() {
    return 1L << mLengthBits;
  }

  /**
   * The length is always a power of 2.
   * length() == 2 ^ bits()
   * @return log_2 of the length.
   */
  public int bits() {
    return mLengthBits;
  }

  /**
   * Get the memory that will be used by the <code>HashBitVector</code>.
   * @return the memory that will be used by the <code>HashBitVector</code>.
   */
  public long bytes() {
    final long lm = length() + BitVector.MASK;
    final long entries = lm >> BitVector.BITS_PER_ENTRY;
    return entries * BitVector.BYTES_PER_ENTRY;
  }

  @Override
  public void toString(final StringBuilder sb) {
    sb.append("HashBit Handle hash bits=").append(mHashBits).append(" length bits=").append(mLengthBits);
  }

  @Override
  public boolean integrity() {
    Exam.assertTrue(mHashBits >= 0 && mHashBits <= 64);
    Exam.assertTrue(mLengthBits >= 0 && mLengthBits <= 62);
    Exam.assertEquals(length(), 1L << bits());
    return true;
  }

}

