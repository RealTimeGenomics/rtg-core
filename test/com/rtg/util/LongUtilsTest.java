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
