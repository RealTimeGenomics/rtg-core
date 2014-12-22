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
package com.rtg.util.array.packedindex;

import com.rtg.util.MathUtils;
import com.rtg.util.array.AbstractIndex;
import com.rtg.util.array.bitindex.BitIndex.IndexType;
import com.rtg.util.array.longindex.LongIndex;
import com.rtg.util.format.FormatInteger;
import com.rtg.util.integrity.Exam;


/**
 * This implements an array of small unsigned integer values,
 * packing as many values as possible into an array of longs.
 * It uses every bit, so a few values will be split across two longs.
 * This is done so that no divisions or modulo operations are needed,
 * and it also minimizes the space used.
 *
 */
public final class PackedIndex extends AbstractIndex {

  /** The main array, in which all values are stored. */
  private final LongIndex mArray;

  /** Number of bits in a long. */
  private static final int LONG_BITS = 64;

  /** The number of bits you need to shift by to divide by 64. */
  private static final int DIVIDE_BY_64 = 6;

  private static final int MODULO_64 = 63;

  /** <code>mLength</code> field is inherited from LongIndex */

  /** Number of values that can be stored in each bit-field. */
  private final long mRange;

  /** Width of each bit-fields */
  private final int mBits;

  /** Equals 2 to the power of <code>mBits</code> minus 1. */
  private final long mMask;

  /**
   * @param length number of items.
   * @param range number of values to be stored in each item.
   * @param type the type of array to use inside the index
   * @exception NegativeArraySizeException if length less than 0
   */
  public PackedIndex(final long length, final long range, final IndexType type) {
    super(length);
    if (range < 2) {
      throw new IllegalArgumentException("Illegal range value=" + range);
    }
    mRange = range;
    mBits = MathUtils.ceilPowerOf2Bits(range - 1);
    mMask = (1L << mBits) - 1;
    final long llen = length * mBits / LONG_BITS + 1;
    // always use LongChunks, otherwise the JIT gets confused.
    mArray = new com.rtg.util.array.longindex.LongChunks(llen);
//    mArray = type == IndexType.DEFAULT ? LongCreate.createIndex(llen)
//        : type == IndexType.CHUNKED ? new com.rtg.util.array.longindex.LongChunks(llen)
//        : new com.rtg.util.array.longindex.Array(llen);
  }

  /**
   * @return the number of bytes consumed.
   */
  @Override
  public long bytes() {
    return mArray.bytes();
  }

  /**
   * Get the number of bits that each value is packed into.
   * @return an integer from 1 to 31.
   */
  public int getBits() {
    return mBits;
  }

  @Override
  public boolean integrity() {
    Exam.assertTrue(mRange <= MathUtils.round(Math.pow(2, mBits)));
    assert mLength >= 0;
    return true;
  }

  @Override
  public long get(final long index) {
    final long bitPos = index * mBits;
    final int  shift = (int) bitPos & MODULO_64;
    final long firstLong = bitPos >> DIVIDE_BY_64;
    final long secondLong = (bitPos + mBits) >> DIVIDE_BY_64;
    final long long1 = mArray.get(firstLong);
    if (firstLong == secondLong) {
      // only need to read one long
      return (long1 >>> shift) & mMask;
    } else {
      final long long2 = mArray.get(secondLong); // might be the same long, but we don't care.
      // this is tricky: a good example of how it works is mBits=7 and shift=60.
      return ((long2 << (64 - shift)) | (long1 >>> shift)) & mMask;
    }
  }

  /**
   * Set the <code>long</code> at the specified index
   *
   * @param index the index
   * @param value the value
   * @throws UnsupportedOperationException if the underlying type
   * is not a <code>long</code>.
   */
  @Override
  public void set(final long index, final long value) {
    assert 0 <= value && value < mRange;
    final long bitPos = index * mBits;
    final int  shift = (int) bitPos & MODULO_64;
    final long firstLong = bitPos >> DIVIDE_BY_64;
    final long secondLong = (bitPos + mBits) >> DIVIDE_BY_64;
    final long long1 = mArray.get(firstLong);
    if (firstLong == secondLong) {
      // it all fits into one long
      final long newLong1 = (long1 & ~(mMask << shift)) | (value << shift);
      mArray.set(firstLong, newLong1);
    } else {
      // it spans two different longs
      final long long2 = mArray.get(secondLong);
      final int  shift2 = 64 - shift;
      final long newLong1 = (long1 & ~(mMask << shift)) | (value << shift);
      final long newLong2 = (long2 & ~(mMask >> shift2)) | (value >> shift2);
      mArray.set(firstLong, newLong1);
      mArray.set(secondLong, newLong2);
    }
  }

  @Override
  protected FormatInteger formatValue() {
    if (mBits <= 16) {
      return new FormatInteger(5);
    }
    if (mBits <= 32) {
      return new FormatInteger(10);
    }
    return new FormatInteger(20);
  }

  @Override
  public boolean safeFromWordTearing() {
    return false;
  }
}

