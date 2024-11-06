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

import junit.framework.TestCase;

/**
 */
public class LongUtilsTest extends TestCase {

  public void testLongMask() {
    for (int i = 1; i < Long.SIZE; ++i) {
      final long mask = LongUtils.longMask(i);
      assertTrue(i + "", 0 == (mask & (1L << i)));
    }
    assertEquals(0L, LongUtils.longMask(0));
    assertEquals(1L, LongUtils.longMask(1));
    assertEquals(3L, LongUtils.longMask(2));
    assertEquals(-1L, LongUtils.longMask(Long.SIZE));
    assertEquals(Long.MAX_VALUE, LongUtils.longMask(Long.SIZE - 1));
  }

  public void testIsLessThanUnsigned() {
    checkIsLessThanUnsigned(new long[] {0L, 1L, 10L, Long.MAX_VALUE, Long.MIN_VALUE, -10L, -1L});
  }

  private void checkIsLessThanUnsigned(long[] ls) {
    for (int i = 0; i < ls.length - 1; ++i) {
      for (int j = i + 1; j < ls.length; ++j) {
        assertTrue(LongUtils.isLessThanUnsigned(ls[i], ls[j]));
        assertFalse(LongUtils.isLessThanUnsigned(ls[j], ls[i]));
      }
      assertFalse(LongUtils.isLessThanUnsigned(ls[i], ls[i]));
    }
  }

  public void testHash() {
    assertEquals(0, LongUtils.hashCode(0));
    assertEquals(1, LongUtils.hashCode(1));
    assertEquals(Integer.MAX_VALUE, LongUtils.hashCode(Integer.MAX_VALUE));
    assertEquals(0, LongUtils.hashCode(-1));
  }

}
