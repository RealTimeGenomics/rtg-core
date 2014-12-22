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
package com.rtg.util.array.zeroindex;


/**
 * Contains the only public ways of constructing a <code>BitIndex</code>.
 */
public final class ZeroCreate {
  private ZeroCreate() { // private so cannot create an instance of this utility class
  }

  /**
   * Create a new ZeroIndex of the specified length with all values 0.
   * @param length number of entries in the index.
   * @return a ZeroIndex.
   * @exception NegativeArraySizeException if length negative.
   * @exception IllegalArgumentException if range is less than 2 or too big.
   */
  public static ZeroIndex createIndex(final long length) throws NegativeArraySizeException, IllegalArgumentException {
    return new ZeroIndex(length, 0);
  }
}
