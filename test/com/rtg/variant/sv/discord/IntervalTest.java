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

package com.rtg.variant.sv.discord;

import com.rtg.util.TestUtils;

import junit.framework.TestCase;

/**
 */
public class IntervalTest extends TestCase {

  public void test1() {
    final Interval ii = new Interval(1, 10);
    ii.integrity();
    assertEquals(1, ii.getA());
    assertEquals(10, ii.getB());
    assertEquals(System.identityHashCode(ii), ii.hashCode());

    final Interval ni = ii.negative();
    ni.integrity();
    assertEquals(-1, ni.getA());
    assertEquals(-10, ni.getB());
  }

  public void test2() {
    final Interval ii = new Interval(-1, -10);
    ii.integrity();
    assertEquals(-1, ii.getA());
    assertEquals(-10, ii.getB());

    final Interval ni = ii.negative();
    ni.integrity();
    assertEquals(1, ni.getA());
    assertEquals(10, ni.getB());
  }

  public void testEquals() {
    TestUtils.equalsTest(new Interval[] {
        new Interval(-1, 1),
        new Interval(0, 1),
        new Interval(-1, 0)
    });
    assertEquals(new Interval(-1, 1), new Interval(-1, 1));
  }
}
