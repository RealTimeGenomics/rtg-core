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

package com.rtg.variant.eval;

import junit.framework.TestCase;

/**
 * Test the corresponding class
 */
public class OrientedVariantTest extends TestCase {

  public void test() {
    final Variant v = new MockVariant(0, 0, new byte[] {0, 1, 2}, null);
    final OrientedVariant ov = new OrientedVariant(v, true);
    assertEquals("0:0 NAC +", ov.toString());
    assertTrue(ov.equals(ov));
    assertEquals(-1, ov.getStart());
    assertEquals(-1, ov.getEnd());
    assertEquals(3, ov.nt(true).length);
    assertEquals(3, ov.ntAlleleA().length);
    assertNull(ov.ntAlleleB());
    final OrientedVariant ov2 = new OrientedVariant(new MockVariant(1, 2, new byte[] {1, 1, 1}, new byte[] {2, 2, 2}), false);
    assertEquals("1:2 AAA:CCC -", ov2.toString());
    assertFalse(ov.equals(ov2));
    assertFalse(ov.hashCode() == ov2.hashCode());
    assertTrue(ov.isAlleleA());
    assertFalse(ov2.isAlleleA());
    assertEquals(0, ov.compareTo(ov));
    assertEquals(-1, ov.compareTo(ov2));
    assertEquals(1, ov2.compareTo(ov));
    assertFalse(ov.equals("not an OrientedVariant"));
    assertFalse(ov.equals(null));

    final Variant v2 = new MockVariant(0, 0, new byte[] {0, 1, 2}, new byte[] {1});
    OrientedVariant ov3 = new OrientedVariant(v2, true);
    OrientedVariant ov4 = new OrientedVariant(v2, false);
    assertEquals(1, ov3.compareTo(ov4));
    assertEquals(-1, ov4.compareTo(ov3));
  }
}
