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

