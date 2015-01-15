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
package com.rtg.util;

/**
 * Utilities for manipulating longs.
 */
public final class LongUtils {

  private LongUtils() { }

  /**
   * Create mask for low order b bits.
   * The tricky bit is getting this right for both <code>b == 64</code> and <code>b == 0</code>;
   * @param b the number of bits.
   * @return a mask which includes <code> 0 ... b-1</code> inclusive.
   */
  public static long longMask(final int b) {
    assert Long.SIZE >= b && b >= 0;
    if (b == 0) {
      return 0;
    }
    return ~0L >>> (Long.SIZE - b);
  }

  /**
   * Check if <code>n1</code> is less than <code>n2</code> when representing 64 bit unsigned longs.
   * @param n1 first long to be compared.
   * @param n2 second long to be compared.
   * @return true iff <code>n1</code> &lt; <code>n2</code>.
   */
  public static boolean isLessThanUnsigned(long n1, long n2) {
    // Note Java 8 offers Long.compareUnsigned that could replace this
    return (n1 < n2) ^ ((n1 < 0) != (n2 < 0));
  }

  /**
   * Return a hashcode for a long.
   * @param v value
   * @return hash
   */
  public static int hashCode(final long v) {
    return (int) (v ^ (v >>> 32));
  }
}
