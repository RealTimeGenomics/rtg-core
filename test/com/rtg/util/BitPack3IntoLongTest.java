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
public class BitPack3IntoLongTest extends TestCase {

  public void test() {
    try {
      new BitPack3IntoLong(32, 32, 1);
      fail();
    } catch (final IllegalArgumentException e) {
    }
    try {
      new BitPack3IntoLong(1, 0, 1);
      fail();
    } catch (final IllegalArgumentException e) {
    }
    try {
      new BitPack3IntoLong(0, 1, 1);
      fail();
    } catch (final IllegalArgumentException e) {
    }
    try {
      new BitPack3IntoLong(1, 1, 0);
      fail();
    } catch (final IllegalArgumentException e) {
    }
    final BitPack3IntoLong bit = new BitPack3IntoLong(15, 16, 19);
    long value = bit.packValues(53, 65535, 524287);
    assertEquals(53, bit.getField(0, value));
    assertEquals(65535, bit.getField(1, value));
    assertEquals(524287, bit.getField(2, value));
    value = bit.packValues(12, 12, 12);
    assertEquals(12, bit.getField(0, value));
    assertEquals(12, bit.getField(1, value));
    assertEquals(12, bit.getField(2, value));
    value = bit.packValues(205, 1105, 405);
    assertEquals(205, bit.getField(0, value));
    assertEquals(1105, bit.getField(1, value));
    assertEquals(405, bit.getField(2, value));
  }
}
