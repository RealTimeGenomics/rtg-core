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
package com.rtg.util.intervals;

import junit.framework.TestCase;

/**
 */
public class StatusIntervalTest extends TestCase {

  private static final byte FOO = 1;
  private static final byte BAR = 2;

  public void test() {
    try {
      new StatusInterval(5, 5);
      fail();
    } catch (IllegalArgumentException e) {
      // ok
    }
    final StatusInterval i = new StatusInterval(5, 10);
    for (int k = 5; k < 10; k++) {
      assertFalse(i.contains(k));
    }
    try {
      i.contains(4);
    } catch (ArrayIndexOutOfBoundsException e) {
      // Expected
    }
    try {
      i.contains(10);
    } catch (ArrayIndexOutOfBoundsException e) {
      // Expected
    }
    i.add(2, 3, FOO);
    i.add(7, 8, FOO);
    for (int k = 5; k < 10; k++) {
      assertTrue((k == 7) ^ !i.contains(k));
    }
    i.add(0, 100, FOO);
    for (int k = 5; k < 10; k++) {
      assertTrue(i.contains(k));
      assertTrue(i.get(k) == FOO);
    }
    i.add(9, 10, BAR);
    assertTrue(i.contains(9));
    assertTrue(i.get(9) == BAR);

  }
}
