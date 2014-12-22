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

/**
 * Contains the only public ways of constructing a ShortIndex.
 */
public final class ByteCreate {
  private ByteCreate() { // private so cant create an instance of this utility class
  }

  /**
   * Create a new ByteIndex of the specified length.
   * Chooses an appropriate implementation depending on the length.
   * @param length number of entries in the ShortIndex.
   * @return a ShortIndex.
   * @exception NegativeArraySizeException if length negative.
   */
  public static ByteIndex createIndex(final long length) {
    if (length < 0) {
      throw new NegativeArraySizeException("Negative length=" + length);
    }
    if (length <= ByteIndex.MAX_LENGTH) {
      return new ByteArray(length);
    } else {
      return new ByteChunks(length);
    }
  }
}

