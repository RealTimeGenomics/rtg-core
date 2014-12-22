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
package com.rtg.util.array.bitindex;


/**
 * Contains the only public ways of constructing a <code>BitIndex</code>.
 */
public final class BitCreate {
  private BitCreate() { // private so cannot create an instance of this utility class
  }

  /**
   * Create a new PackedIndex of the specified length,
   * which can store unsigned values <code>0 .. range-1</code>.
   * @param length number of entries in the index.
   * @param bits the number of bits required to store each value.
   * @return a PackedIndex.
   * @exception NegativeArraySizeException if length negative.
   * @exception IllegalArgumentException if range is less than 2 or too big.
   */
  public static BitIndex createIndex(final long length, final int bits) throws NegativeArraySizeException, IllegalArgumentException {
    if (bits < 1 || bits > 64) {
      throw new IllegalArgumentException("Illegal bits value=" + bits);
    }
    return new BitIndex(length, bits);
  }
}

