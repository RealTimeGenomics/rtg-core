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
public class BitPackHelperLongTest extends TestCase {

  public BitPackHelperLongTest(String testName) {
    super(testName);
  }

  public void test() {
    BitPackHelperLong bit;
    try {
      new BitPackHelperLong(32, 16, 16, 1);
      fail();
    } catch (final IllegalArgumentException e) {

    }
    bit = new BitPackHelperLong(14, 15, 16, 19);
    long value = bit.packValues(new long[] {34, 53, 65535, 524287});
    assertEquals(34, bit.getField(0, value));
    assertEquals(53, bit.getField(1, value));
    assertEquals(65535, bit.getField(2, value));
    assertEquals(524287, bit.getField(3, value));
    value = bit.packValues(new long[] {12, 12, 12, 12});
    assertEquals(12, bit.getField(0, value));
    assertEquals(12, bit.getField(1, value));
    assertEquals(12, bit.getField(2, value));
    assertEquals(12, bit.getField(3, value));
    value = bit.packValues(new long[] {205, 1105, 405, 905});
    assertEquals(205, bit.getField(0, value));
    assertEquals(1105, bit.getField(1, value));
    assertEquals(405, bit.getField(2, value));
    assertEquals(905, bit.getField(3, value));
  }
}
