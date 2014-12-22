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
 * Tests for <code>Resources</code>.
 *
 */
public class ReverseDoubleTest extends TestCase {

  public void test() {
    final ReverseDouble rd = new ReverseDouble(1.0);
    assertEquals(1.0, rd.doubleValue());
    assertEquals(Double.toString(1.0), rd.toString());
    assertFalse(rd.equals(null));
    assertFalse(rd.equals(new Object()));
  }

  public void testHash() {
    final double dbl = Double.longBitsToDouble(1L << 32);
    final ReverseDouble rd = new ReverseDouble(dbl);
    assertEquals(1, rd.hashCode());
  }

  public final void testOrder() {
    TestUtils.testOrder(new ReverseDouble[] {
        new ReverseDouble(Double.POSITIVE_INFINITY),
        new ReverseDouble(1000.0),
        new ReverseDouble(1.0),
        new ReverseDouble(0.0),
        new ReverseDouble(-1.0),
        new ReverseDouble(-1000.0),
        new ReverseDouble(Double.NEGATIVE_INFINITY),
    }, true);
  }
}

