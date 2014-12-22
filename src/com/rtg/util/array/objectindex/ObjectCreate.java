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
package com.rtg.util.array.objectindex;

/**
 * Contains the only public ways of constructing an IntIndex.
 */
public final class ObjectCreate {
  private ObjectCreate() { // private so cant create an instance of this utility class
  }

  /**
   * Create a new IntIndex of the specified length.
   * Chooses an appropriate implementation depending on the length.
   * @param length number of entries in the IntIndex.
   * @return an IntIndex.
   * @param <A> type stored in <code>ObjectIndex</code>
   * @exception NegativeArraySizeException if length negative.
   */
  public static <A> ObjectIndex<A> createIndex(final long length) {
    if (length < 0) {
      throw new NegativeArraySizeException("Negative length=" + length);
    }
    if (length <= ObjectIndex.MAX_LENGTH) {
      return new ObjectArray<>(length);
    } else {
      return new ObjectChunks<>(length);
    }
  }
}

