/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
 */

package com.rtg.index;

import java.util.Arrays;

import com.rtg.AbstractTest;

/**
 * Test
 */
public class UnionRepeatFrequencyFilterMethodTest extends AbstractTest {

  public void test() {
    final IndexFilterMethod m = new UnionRepeatFrequencyFilterMethod(
      new FixedRepeatFrequencyFilterMethod(50),
      new BlacklistFilterMethod(Arrays.asList(0b1111L, 0b1010L, 0b1100L), 4, 1));
    m.initialize(null);
    // Fail blacklist, (freq OK)
    assertFalse(m.keepHash(0b1111L, 10));
    assertFalse(m.keepHash(0b1100L, 10));
    assertFalse(m.keepHash(0b1010L, 10));
    // Fail frequency, (blacklist OK)
    assertFalse(m.keepHash(0, 80));
    assertFalse(m.keepHash(0, 1000));
    assertFalse(m.keepHash(0, 10000));
    assertFalse(m.keepHash(0b1011L, 60));
    // Pass both
    assertTrue(m.keepHash(0b1110L, 10));
    assertTrue(m.keepHash(0b1101L, 10));
    assertTrue(m.keepHash(0b1011L, 10));
    assertTrue(m.keepHash(0, 5));
    assertTrue(m.keepHash(0, 25));
    assertTrue(m.keepHash(0, 50));
  }
}