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

package com.rtg.variant.sv;

import junit.framework.TestCase;

/**
 */
public class DistributionStepTest extends TestCase {

  public void test() {
    final Distribution da = new DistributionStep(-5, 5, 2, 2.0, 0.1);
    da.globalIntegrity();
    assertEquals(-5, da.lo());
    assertEquals(5, da.hi());
    assertEquals(2.0, da.get(-5));
    assertEquals(2.0, da.get(-1));
    assertEquals(2.0, da.get(0));
    assertEquals(2.0, da.get(1));
    assertEquals(2.0, da.get(2));
    assertEquals(0.1, da.get(3));
    assertEquals(0.1, da.get(4));

    assertEquals("Step:radius=5 offset=2 rate1=2.0000 rate2=0.1000", da.toString());

    try {
      da.get(-6);
      fail();
    } catch (final IndexOutOfBoundsException e) {
      assertEquals("index=-6 lo=-5 hi=5", e.getMessage());
    }
    try {
      da.get(5);
      fail();
    } catch (final IndexOutOfBoundsException e) {
      assertEquals("index=5 lo=-5 hi=5", e.getMessage());
    }
  }

}
