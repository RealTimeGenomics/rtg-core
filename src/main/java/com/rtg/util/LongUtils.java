/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
    return b == 0 ? 0 : ~0L >>> (Long.SIZE - b);
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
