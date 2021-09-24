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
package com.rtg.variant.cnv.region;

import com.rtg.util.TestUtils;

import junit.framework.TestCase;

/**
 */
public class DummyCnvRegionTest extends TestCase {

  public void test() {
    final AbstractCnvRegion r = new SimpleCnvRegion(0, 1);
    assertTrue(r.contains(0));
    assertFalse(r.contains(-1));
    assertFalse(r.contains(1));
    assertEquals(0, r.getStart());
    assertEquals(1, r.getEnd());
  }

  public void testBad() {
    try {
      new SimpleCnvRegion(1, 0);
      fail();
    } catch (final IllegalArgumentException e) {
      //expected
    }
  }

  public void testComparable() {
    final AbstractCnvRegion r0 = new SimpleCnvRegion(0, 0);
    final AbstractCnvRegion r1 = new SimpleCnvRegion(1, 1);
    final AbstractCnvRegion r2 = new SimpleCnvRegion(5, 7);
    final AbstractCnvRegion r3 = new SimpleCnvRegion(5, 8);
    final AbstractCnvRegion r4 = new SimpleCnvRegion(6, 23);
    final AbstractCnvRegion r5 = new SimpleCnvRegion(20, Integer.MAX_VALUE);
    TestUtils.testOrder(new AbstractCnvRegion[][] {{r0, r0}, {r1, r1}, {r2, r2}, {r3, r3}, {r4, r4}, {r5, r5}}, false);
  }
}
