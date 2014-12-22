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
package com.rtg.util.array.byteindex;

import com.rtg.util.array.AbstractIndex;
import com.rtg.util.format.FormatInteger;

/**
 * Common code used in implementing all the byte index variants. Holds some handy
 * constants as well as the length of the index.
 *
 */
public abstract class ByteIndex extends AbstractIndex {

  /** Number of bits in a short. */
  static final int BYTE_BITS = 8;

  /** The low order bits of a long corresponding to a short. */
  static final long BYTE_MASK = (1L << BYTE_BITS) - 1L;

  /** The bits above those used by a short. */
  static final long HIGH_MASK = ~BYTE_MASK;

  /** The bits from the signed bit for a short up. */
  static final long HIGH_SIGNED_MASK = ~((1L << (BYTE_BITS - 1)) - 1L);

  /**
   * Maximum number of bits that can be used when allocating a byte array.
   */
  protected static final int MAX_BITS = 30;

  /**
   * Length of largest allocatable short array.
   */
  static final long MAX_LENGTH = 1L << MAX_BITS;

  /**
   * Information used in creating "chunks" in some of the implementations. Be
   * wary of changing CHUNK_BITS.
   */
  protected static final int CHUNK_BITS = 30;

  /**
   * Number of bytes in a short.
   */
  protected static final int BYTE_SIZE = 1;

  /**
   * @param length of the array.
   * @exception NegativeArraySizeException if length less than 0
   */
  public ByteIndex(final long length) {
    super(length);
  }

  /**
   * Swap the values at the two specified locations.
   *
   * @param index1 the first index to be swapped
   * @param index2 the second index to be swapped
   */
  @Override
  public void swap(final long index1, final long index2) {
    // Default implementation - can often be made faster in particular
    // implementations
    final byte temp = getByte(index1);
    setByte(index1, getByte(index2));
    setByte(index2, temp);
  }

  /**
   * @return the number of bytes consumed.
   */
  @Override
  public long bytes() {
    return BYTE_SIZE * mLength;
  }

  static final FormatInteger FORMAT_VALUE = new FormatInteger(5);

  @Override
  protected FormatInteger formatValue() {
    return FORMAT_VALUE;
  }

  @Override
  public final void set(final long index, final long value) {
    //High order bits must be zero
    assert (value & HIGH_MASK) == 0L : value;
    setByte(index, (byte) value);
  }

  @Override
  public final long get(final long index) {
    return getByte(index) & BYTE_MASK; //clear any propogated sign bits
  }

  /**
   * Get the <code>byte</code> at the specified index
   *
   * @param index the index
   * @return byte value
   * @throws UnsupportedOperationException if the underlying type
   * is not a <code>long</code>.
   */
  public abstract byte getByte(final long index);

  /**
   * Set the <code>byte</code> at the specified index
   *
   * @param index the index
   * @param value the value
   * @throws UnsupportedOperationException if they underlying type
   * is not a <code>long</code>.
   */
  public abstract void setByte(final long index, final byte value);

}

