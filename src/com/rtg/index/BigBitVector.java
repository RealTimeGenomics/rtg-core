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

import com.rtg.util.array.intindex.IntChunks;
import com.rtg.util.integrity.Exam;


/**
 * A vector of bits.
 * Uses underlying big arrays so that can be longer than 2^34 bits.
 */
public final class BigBitVector extends AbstractBitVector {

  private final IntChunks mArray;

  /**
   * Creates a new <code>BitVector</code> with the specified length.
   *
   * @param length the number of bits.
   */
  public BigBitVector(final long length) {
    super(length);
    final long allocLength = (length + MASK) >> BITS_PER_ENTRY;
    mArray = new IntChunks(allocLength);
  }

  @Override
  public boolean get(final long index) {
    if (index >= mLength || index < 0) {
      throw new ArrayIndexOutOfBoundsException(index + ":" + mLength);
    }
    final long x = index >> BITS_PER_ENTRY;
    final int v = mArray.getInt(x);
    final int s = 1 << (int) (index & MASK);
    return (v & s) != 0;
  }

  @Override
  public void set(final long index) {
    if (index >= mLength || index < 0) {
      throw new ArrayIndexOutOfBoundsException(index + ":" + mLength);
    }
    final long x = index >> BITS_PER_ENTRY;
    final int v = mArray.getInt(x) | 1 << (int) (index & MASK);
    mArray.setInt(x, v);
  }

  @Override
  public void reset(final long index) {
    if (index >= mLength || index < 0) {
      throw new ArrayIndexOutOfBoundsException(index + ":" + mLength);
    }
    final long x = index >> BITS_PER_ENTRY;
    final int v = mArray.getInt(x) ^ (1 << (int) (index & MASK));
    mArray.setInt(x, v);
  }


  @Override
  public long bytes() {
    return (mArray.length() << BITS_PER_ENTRY) >> BITS_PER_BYTE;
  }


  /**
   * Gets the length of the bit vector.
   *
   * @return number of bits.
   */
  @Override
  public long length() {
    return mLength;
  }

  @Override
  public boolean integrity() {
    if (mArray == null) {
      Exam.assertTrue(false);
    } else {
      Exam.assertTrue(mArray.length() << BITS_PER_ENTRY >= mLength);
    }
    return true;
  }

}


