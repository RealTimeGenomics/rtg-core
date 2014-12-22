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
 *
 */
public class BitPack3IntoIntTest extends TestCase {

  public BitPack3IntoIntTest(String testName) {
    super(testName);
  }

  public void test() {
    BitPack3IntoInt bit;
    try {
      new BitPack3IntoInt(32, 1, 0);
      fail();
    } catch (final IllegalArgumentException e) {

    }
    bit = new BitPack3IntoInt(13, 8, 11);
    int value = bit.packValues(34, 90, 1000);
    assertEquals(34, bit.getField(0, value));
    assertEquals(90, bit.getField(1, value));
    assertEquals(1000, bit.getField(2, value));
    value = bit.packValues(12, 12, 12);
    assertEquals(12, bit.getField(0, value));
    assertEquals(12, bit.getField(1, value));
    assertEquals(12, bit.getField(2, value));
    value = bit.packValues(25, 11, 505);
    assertEquals(25, bit.getField(0, value));
    assertEquals(11, bit.getField(1, value));
    assertEquals(505, bit.getField(2, value));
  }

}
