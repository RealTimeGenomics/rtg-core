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
public class BitPackHelperIntTest extends TestCase {

  public BitPackHelperIntTest(String testName) {
    super(testName);
  }

  public void test() {
    BitPackHelperInt bit;
    try {
      new BitPackHelperInt(32, 1);
      fail();
    } catch (final IllegalArgumentException e) {

    }
    bit = new BitPackHelperInt(6, 7, 8, 11);
    int value = bit.packValues(new int[] {34, 90, 255, 1000});
    assertEquals(34, bit.getField(0, value));
    assertEquals(90, bit.getField(1, value));
    assertEquals(255, bit.getField(2, value));
    assertEquals(1000, bit.getField(3, value));
    value = bit.packValues(new int[] {12, 12, 12, 12});
    assertEquals(12, bit.getField(0, value));
    assertEquals(12, bit.getField(1, value));
    assertEquals(12, bit.getField(2, value));
    assertEquals(12, bit.getField(3, value));
    value = bit.packValues(new int[] {25, 11, 94, 505});
    assertEquals(25, bit.getField(0, value));
    assertEquals(11, bit.getField(1, value));
    assertEquals(94, bit.getField(2, value));
    assertEquals(505, bit.getField(3, value));
  }
}
