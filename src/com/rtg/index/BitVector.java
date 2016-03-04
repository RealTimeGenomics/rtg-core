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
/**
 * A vector of bits.
 */
public final class BitVector extends AbstractBitVector {

  private final int[] mArray;

  /**
   * Creates a new <code>BitVector</code> with the specified length.
   *
   * @param length the number of bits.
   */
  public BitVector(final long length) {
    super(length);
    final long allocLength = (length + MASK) >> BITS_PER_ENTRY;
    if (allocLength >= (Integer.MAX_VALUE >> 2)) {
      throw new RuntimeException("BitVector is too long " + length);
    }
    mArray = new int[(int) allocLength];
  }

  private void checkBounds(final long index) {
    if (index >= mLength || index < 0) {
      throw new ArrayIndexOutOfBoundsException(index + ":" + mLength);
    }
  }

  @Override
  public boolean get(final long index) {
    checkBounds(index);
    final long x = index >> BITS_PER_ENTRY;
    final int v = mArray[(int) x];
    final int s = 1 << (int) (index & MASK);
    return (v & s) != 0;
  }

  @Override
  public void set(final long index) {
    checkBounds(index);
    final long x = index >> BITS_PER_ENTRY;
    mArray[(int) x] |= 1 << (int) (index & MASK);
  }

  @Override
  public void reset(final long index) {
    checkBounds(index);
    final long x = index >> BITS_PER_ENTRY;
    mArray[(int) x] &= ~(1 << (int) (index & MASK));
  }

  @Override
  public long bytes() {
    return (((long) mArray.length) << BITS_PER_ENTRY) >> BITS_PER_BYTE;
  }

  @Override
  public long length() {
    return mLength;
  }
}


