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
package com.rtg.index.hash.ngs.general;


import junit.framework.TestCase;

/**
 */
public class UtilTest extends TestCase {

  /**
   * Test method for {@link com.rtg.index.hash.ngs.general.Util#binomial(int, int)}.
   */
  public final void testBinomial() {
    assertEquals(1, Util.binomial(0, 0));
    assertEquals(0, Util.binomial(0, -1));
    assertEquals(0, Util.binomial(0, 1));
    assertEquals(0, Util.binomial(1, -1));
    assertEquals(1, Util.binomial(1, 0));
    assertEquals(1, Util.binomial(1, 1));
    assertEquals(0, Util.binomial(1, 2));
    assertEquals(0, Util.binomial(2, -1));
    assertEquals(1, Util.binomial(2, 0));
    assertEquals(2, Util.binomial(2, 1));
    assertEquals(1, Util.binomial(2, 2));
    assertEquals(0, Util.binomial(2, 3));
    assertEquals(10, Util.binomial(5, 3));
    assertEquals((63 * 62) / 2, Util.binomial(63, 2));
  }

  public void testBad() {
    try {
      Util.binomial(65, 3);
      fail();
    } catch (final ArrayIndexOutOfBoundsException e) {
      //expected
    }
  }

}
