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
package com.rtg.util.array.longindex;

import java.io.IOException;
import java.io.ObjectOutputStream;

import com.rtg.util.array.AbstractIndex;
import com.rtg.util.format.FormatInteger;

/**
 * Common code used in implementing all the long index variants. Holds some handy
 * constants as well as the length of the index.
 *
 */
public abstract class LongIndex extends AbstractIndex {

  /**
   * Information used in creating "chunks" in some of the implementations. Be
   * wary of changing CHUNK_BITS.
   */
  protected static final int CHUNK_BITS = 27;

  /**
   * Number of bytes in a long.
   */
  protected static final int LONG_SIZE = 8;

  /**
   * @param length of the array.
   * @exception NegativeArraySizeException if length less than 0
   */
  protected LongIndex(final long length) {
    super(length);
  }

  /**
   * @return the number of bytes consumed.
   */
  @Override
  public long bytes() {
    return LONG_SIZE * mLength;
  }

  static final FormatInteger FORMAT_VALUE = new FormatInteger(20);

  @Override
  protected FormatInteger formatValue() {
    return FORMAT_VALUE;
  }

  /**
   * Save this index such that it can be loaded again from {@link LongCreate#loadIndex(java.io.ObjectInputStream)}
   * @param dos steam to save to
   * @throws IOException if an IO error occurs
   */
  public abstract void save(ObjectOutputStream dos) throws IOException;
}
