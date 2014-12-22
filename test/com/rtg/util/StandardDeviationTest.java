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

import java.util.Arrays;

import junit.framework.TestCase;

/**
 */
public class StandardDeviationTest extends TestCase {
  public void testStandardDeviation() {
    StandardDeviation dev = new StandardDeviation();
    assertEquals(0.0, dev.standardDeviation());
    assertEquals(0.0, dev.mean());
    dev.addSample(3);
    dev.addSample(4);
    assertEquals(0.5, dev.standardDeviation());
    assertEquals(3.5, dev.mean());
    assertEquals("2\t3.5\t0.5\t7.0\t25.0", dev.toString());
  }
  public void test2() {
    final StandardDeviation dev = new StandardDeviation();
    for (int i : new int[] {2, 4, 4, 4, 5, 5, 7, 9}) {
      dev.addSample(i);
    }
    assertEquals(2.0, dev.standardDeviation(), 0.1);
    assertEquals(5.0, dev.mean(), 0.1);
  }

  public void testCombine() {
    StandardDeviation first = new StandardDeviation();
    StandardDeviation second = new StandardDeviation();
    for (int i : new int[] {2, 4, 4, 4}) {
      first.addSample(i);
    }
    for (int i : new int[] {5, 5, 7, 9}) {
      second.addSample(i);
    }
    StandardDeviation dev = StandardDeviation.combine(Arrays.asList(first, new StandardDeviation(), second));
    assertEquals(2.0, dev.standardDeviation(), 0.1);
    assertEquals(5.0, dev.mean(), 0.1);
  }
}
