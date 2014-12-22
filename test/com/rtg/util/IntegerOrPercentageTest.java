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
 * Test Integer Or Percentage
 *
 */
public class IntegerOrPercentageTest extends TestCase {
  public void testConstructor() {
    IntegerOrPercentage p = new IntegerOrPercentage("10%");
    assertTrue(p.isPercentage());
    assertEquals(10, p.getValue(100));
    assertEquals(200, p.getValue(2000));
    assertEquals("10%", p.toString());
    IntegerOrPercentage q = new IntegerOrPercentage("10");
    assertFalse(q.isPercentage());
    assertEquals(10, q.getValue(100));
    assertEquals(10, q.getValue(2000));
    assertEquals("10", q.toString());

    assertTrue(p.compareTo(q) < 0);
    assertTrue(q.compareTo(p) > 0);
    assertTrue(p.compareTo(new IntegerOrPercentage("15%")) < 0);
    assertEquals(0, q.compareTo(q));
    assertEquals(0, q.compareTo(q));

    IntegerOrPercentage r = new IntegerOrPercentage(15);
    assertFalse(r.isPercentage());
    assertEquals(15, r.getValue(100));
    assertEquals(15, r.getValue(2000));

    assertEquals(p, new IntegerOrPercentage("10%"));
    assertEquals(p, IntegerOrPercentage.valueOf("10%"));
    assertEquals(r, IntegerOrPercentage.valueOf(15));

    assertFalse(r.equals(p));
    assertFalse(r.equals("Monkey"));

    assertEquals(10, p.getRawValue());

    // Shuts jumble up doesn't really do much else
    assertEquals(7038, p.hashCode());
    assertEquals(7037, q.hashCode());
  }
}
