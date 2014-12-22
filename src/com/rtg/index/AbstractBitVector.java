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

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.util.StringUtils;
import com.rtg.util.integrity.IntegralAbstract;

/**
 * A vector of bits.
 *
 */
@TestClass(value = {"com.rtg.index.BitVectorTest"})
public abstract class AbstractBitVector extends IntegralAbstract {

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
  public void toString(final StringBuilder sb) {
    sb.append("BitVector[").append(mLength).append("]").append(StringUtils.LS);
    for (long l = 0; l < mLength; l += 100) {
      sb.append("[").append(l).append("]\t");
      final long min = Math.min(l + 100, mLength);
      for (long m = l; m < min; m++) {
        if (m % 10 == 0) {
          sb.append(" ");
        }
        sb.append(get(m) ? "1" : "0");
      }
      sb.append(StringUtils.LS);
    }
    sb.append(StringUtils.LS);
  }
}
