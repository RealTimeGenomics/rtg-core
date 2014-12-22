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

import com.rtg.util.array.bitindex.BitIndex;


/**
 * Contains the only public ways of constructing a PackedIndex.
 */
public final class PackedCreate {
  private PackedCreate() { // private so cannot create an instance of this utility class
  }

  /**
   * Create a new PackedIndex of the specified length,
   * which can store unsigned values <code>0 .. range-1</code>.
   * @param length number of entries in the index.
   * @param range the values can range from 0 up to <code>range - 1</code>.
   * @return a PackedIndex.
   * @exception NegativeArraySizeException if length negative.
   * @exception IllegalArgumentException if range is less than 2 or too big.
   */
  public static PackedIndex createIndex(long length, long range) throws NegativeArraySizeException, IllegalArgumentException {
    return new PackedIndex(length, range, BitIndex.IndexType.DEFAULT);
  }
}

