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

import static com.rtg.util.StringUtils.LS;

import junit.framework.TestCase;

/**
 */
public class DistributionArrayTest extends TestCase {

  public void test() {
    final Distribution da = new DistributionArray(-2, new double[] {1.0, 2.0, 3.0, 4.0});
    da.globalIntegrity();
    assertEquals(-2, da.lo());
    assertEquals(2, da.hi());

    assertEquals(1.0, da.get(-2));
    assertEquals(2.0, da.get(-1));
    assertEquals(3.0, da.get(0));
    assertEquals(4.0, da.get(1));

    assertEquals("[1.000, 2.000, 3.000, 4.000]", da.toString());
    final String exp = ""
      + "-2 1.0000" + LS
      + "-1 2.0000" + LS
      + "0 3.0000" + LS
      + "1 4.0000" + LS
      ;
    assertEquals(exp, da.dump());

    try {
      da.get(-3);
      fail();
    } catch (final IndexOutOfBoundsException e) {
      assertEquals("index=-3 lo=-2 hi=2", e.getMessage());
    }
    try {
      da.get(2);
      fail();
    } catch (final IndexOutOfBoundsException e) {
      assertEquals("index=2 lo=-2 hi=2", e.getMessage());
    }
  }
}
