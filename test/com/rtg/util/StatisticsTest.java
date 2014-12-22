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
public class StatisticsTest extends TestCase {
  public void testSmall() {
    final Statistics stats = new Statistics();
    assertEquals(0.0, stats.standardDeviation());
    assertEquals(0.0, stats.mean());
    assertEquals(0, stats.count());
    assertEquals(0, stats.min());
    assertEquals(0, stats.max());
    assertEquals(0, stats.sum());

    stats.addSample(3);
    stats.addSample(4);
    assertEquals(0.707, stats.standardDeviation(), 0.001);
    assertEquals(3.5, stats.mean());
    assertEquals(2, stats.count());
    assertEquals(3, stats.min());
    assertEquals(4, stats.max());
    assertEquals(7, stats.sum());

    stats.reset();
    assertEquals(0.0, stats.standardDeviation());
    assertEquals(0.0, stats.mean());
    assertEquals(0, stats.count());
    assertEquals(0, stats.min());
    assertEquals(0, stats.max());
    assertEquals(0, stats.sum());
  }

  public void test2() {
    final Statistics stats = new Statistics();
    for (int i : new int[] {2, 4, 4, 4, 5, 5, 7, 9}) {
      stats.addSample(i);
    }
    assertEquals(2.138, stats.standardDeviation(), 0.001);
    assertEquals(5.0, stats.mean(), 0.1);
    assertEquals(8, stats.count());
    assertEquals(2, stats.min());
    assertEquals(9, stats.max());
    assertEquals(40, stats.sum());
  }

  public void test3() {
    final Statistics stats = new Statistics();
    long sum = 0L;
    for (int i = 100; i >= 0; i--) {
      sum += i;
      stats.addSample(i);
      //assertEquals(2.138, stats.standardDeviation(), 0.001);
      assertEquals(sum / (double) (101 - i), stats.mean(), 0.001);
      assertEquals(101 - i, stats.count());
      assertEquals(i, stats.min());
      assertEquals(100, stats.max());
      assertEquals(sum, stats.sum());
    }
    assertEquals(29.300, stats.standardDeviation(), 0.001);
  }

  public void testToString() {
    final Statistics stats = new Statistics();
    for (int i : new int[] {9, 7, 5, 5, 4, 4, 4, 2}) {
      stats.addSample(i);
    }
    assertEquals(2.138, stats.standardDeviation(), 0.001);
    assertEquals(5.0, stats.mean(), 0.1);
    assertEquals(8, stats.count());
    assertEquals(2, stats.min());
    assertEquals(9, stats.max());
    assertEquals(40, stats.sum());

    final String res = stats.toString();
    for (String s : new String[] {"sum  : 40",
        "min  : 2",
        "max  : 9",
        "ss   : 232",
        "n    : 8",
        "mean : 5.0000",
        "sd   : 2.1381"}) {
      assertTrue(s, res.contains(s));
    }
  }
}
