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
public class BitPack2IntoLongTest extends TestCase {

  public void test() {
    try {
      new BitPack2IntoLong(64, 1);
      fail();
    } catch (final IllegalArgumentException e) {
    }
    try {
      new BitPack2IntoLong(60, 0);
      fail();
    } catch (final IllegalArgumentException e) {
    }
    try {
      new BitPack2IntoLong(0, 60);
      fail();
    } catch (final IllegalArgumentException e) {
    }
    final BitPack2IntoLong bit = new BitPack2IntoLong(13, 8);
    long value = bit.packValues(34, 90);
    assertEquals(34, bit.getField(0, value));
    assertEquals(90, bit.getField(1, value));
    value = bit.packValues(12, 12);
    assertEquals(12, bit.getField(0, value));
    assertEquals(12, bit.getField(1, value));
    value = bit.packValues(25, 11);
    assertEquals(25, bit.getField(0, value));
    assertEquals(11, bit.getField(1, value));
  }
}
