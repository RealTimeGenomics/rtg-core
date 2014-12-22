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
package com.rtg.util.array.shortindex;

import java.io.IOException;
import java.io.ObjectOutputStream;

import com.rtg.util.array.AbstractIndex;
import com.rtg.util.format.FormatInteger;

/**
 * Common code used in implementing all the short index variants. Holds some handy
 * constants as well as the length of the index.
 *
 */
public abstract class ShortIndex extends AbstractIndex {

  /** Number of bits in a short. */
  static final int SHORT_BITS = 16;

  /** The low order bits of a long corresponding to a short. */
  static final long SHORT_MASK = (1L << SHORT_BITS) - 1L;

  /** The bits above those used by a short. */
  static final long HIGH_MASK = ~SHORT_MASK;

  /** The bits from the signed bit for a short up. */
  static final long HIGH_SIGNED_MASK = ~((1L << (SHORT_BITS - 1)) - 1L);

  /**
   * Maximum number of bits that can be used when allocating a short array.
   */
  protected static final int MAX_BITS = 29;

  /**
   * Length of largest allocatable short array.
   */
  static final long MAX_LENGTH = 1L << MAX_BITS;

  /**
   * Information used in creating "chunks" in some of the implementations. Be
   * wary of changing CHUNK_BITS.
   */
  protected static final int CHUNK_BITS = 29;

  /**
   * Number of bytes in a short.
   */
  protected static final int SHORT_SIZE = 2;

  /**
   * @param length of the array.
   * @exception NegativeArraySizeException if length less than 0
   */
  protected ShortIndex(final long length) {
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
    final short temp = getShort(index1);
    setShort(index1, getShort(index2));
    setShort(index2, temp);
  }

  /**
   * @return the number of bytes consumed.
   */
  @Override
  public long bytes() {
    return SHORT_SIZE * mLength;
  }

  static final FormatInteger FORMAT_VALUE = new FormatInteger(5);

  @Override
  protected FormatInteger formatValue() {
    return FORMAT_VALUE;
  }

  @Override
  public final void set(final long index, final long value) {
    //High order bits must be zero
    assert (value & HIGH_MASK) == 0L;
    setShort(index, (short) value);
  }

  @Override
  public final long get(final long index) {
    return getShort(index) & SHORT_MASK; //clear any propogated sign bits
  }

  /**
   * Get the <code>short</code> at the specified index
   *
   * @param index the index
   * @return long value
   * @throws UnsupportedOperationException if the underlying type
   * is not a <code>long</code>.
   */
  public abstract short getShort(final long index);

  /**
   * Set the <code>short</code> at the specified index
   *
   * @param index the index
   * @param value the value
   * @throws UnsupportedOperationException if they underlying type
   * is not a <code>long</code>.
   */
  public abstract void setShort(final long index, final short value);

  /**
   * Save this index such that it can be loaded again from {@link ShortCreate#loadIndex(java.io.ObjectInputStream)}
   * @param dos steam to save to
   * @throws IOException if an IO error occurs
   */
  public abstract void save(ObjectOutputStream dos) throws IOException;
}

