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
    mBitVector.toString(sb);
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

