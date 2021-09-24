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

package com.rtg.calibrate;

import junit.framework.TestCase;

/**
 */
public class CovariateSingleReadGroupTest extends TestCase {

  public void test() {
    final CovariateSingleReadGroup rg = new CovariateSingleReadGroup("test");
    assertEquals("readgroup", rg.name());
    assertEquals(CovariateEnum.READGROUP, rg.getType());
    assertEquals(1, rg.size());
    assertEquals(0, rg.parse("test"));
    try {
      rg.parse("foo");
      fail();
    } catch (final UnsupportedOperationException e) {
      // ok
    }
    try {
      rg.valueOf("foo");
      fail();
    } catch (final UnsupportedOperationException e) {
      // ok
    }
    assertEquals(1, rg.newSize());
    assertEquals(1, rg.size());
    assertEquals(false, rg.sizeChanged());
    assertEquals("test", rg.valueString(0));
    assertEquals("test", rg.valueString(1));
  }
}
