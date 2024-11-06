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
public class HashBitVector extends IntegralAbstract {

  /** Number of bits in a Long. */
  private static final int LONG_BITS = 64;

  /** Number of bits in the hash code. */
  private final int mBits;

  /** Number of bits in the index to the bit vector. */
  private final int mVectorBits;

  private final AbstractBitVector mBitVector;

  /**
   * The number of bits to shift the hash in order to get the index into the bit
   * vector.
   */
  private final int mShift;

  /**
   * @param bits the number of valid bits used in the original hash.
   * @param vectorBits the number of bits to be used to address the vector.
   */
  public HashBitVector(final int bits, final int vectorBits) {
    if (vectorBits < 0 || bits < 0 || bits > LONG_BITS) {
      throw new RuntimeException("Invalid bit parameters bits=" + bits + " vectorBits=" + vectorBits);
    }
    mBits = bits;
    mVectorBits = vectorBits;
    final long len = 1L << mVectorBits;
    mBitVector = AbstractBitVector.createBitVector(len);
    mShift = mBits >= mVectorBits ? mBits - mVectorBits : 0;
    //System.err.println("bits=" + mBits + " vectorBits=" + mVectorBits + " shift=" + mShift);
  }

  /**
   * Takes the high order bits of the hash and sets the corresponding bits.
   * @param hash to be set.
   */
  public void set(final long hash) {
    if (mBits != 0) {
      assert mBits == LONG_BITS || hash >>> mBits == 0;
    final long index = hash >>> mShift;
    mBitVector.set(index);
    }
  }

  /**
   * Check if the bit has been set corresponding to the hash.
   * @param hash to be checked.
   * @return true iff the bit corresponding to the hash has been set.
   */
  public boolean get(final long hash) {
    if (mBits != 0) {
      final long index = hash >>> mShift;
      return mBitVector.get(index);
    } else {
      return false;
    }
  }

  /**
   * Set the corresponding bit.
   * @param index bit to be set.
   */
  public void setDirect(final long index) {
    mBitVector.set(index);
  }

  /**
   * Reset the corresponding bit.
   * @param index bit to be reset.
   */
  public void resetDirect(final long index) {
    mBitVector.reset(index);
  }

  /**
   * Check if the bit has been set.
   * @param index to be checked.
   * @return true iff the bit corresponding to the hash has been set.
   */
  public boolean getDirect(final long index) {
    return mBitVector.get(index);
  }

  @Override
  public void toString(final StringBuilder sb) {
    sb.append("HashBitVector bits=").append(mBits).append(" vectorBits=").append(mVectorBits).append(" shift=").append(mShift).append(com.rtg.util.StringUtils.LS);
    sb.append(mBitVector);
  }

  @Override
  public boolean integrity() {
    Exam.assertTrue(mBits > 0 && mBits <= LONG_BITS);
    Exam.assertTrue(mVectorBits >= 0);
    if (mBitVector == null) {
      Exam.assertTrue(false);
    } else {
      Exam.assertEquals(mBitVector.length(), 1L << mVectorBits);
    }
    Exam.assertTrue(mShift >= 0 && mShift < 64);
    return true;
  }

  /**
   * Gets the number of bytes in the vector.
   *
   * @return number of bytes in bit vector.
   */
  public long bytes() {
    return mBitVector.bytes();
  }


  /**
   * Gets the length of the bit vector.
   *
   * @return number of bits.
   */
  public long length() {
    return mBitVector.length();
  }
}

