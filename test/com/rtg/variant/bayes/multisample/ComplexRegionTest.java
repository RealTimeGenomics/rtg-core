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
package com.rtg.variant.bayes.multisample;

import com.rtg.variant.bayes.multisample.ComplexRegion.RegionType;

import junit.framework.TestCase;

/**
 * Test class
 */
public class ComplexRegionTest extends TestCase {

  public void test() {
    check(ComplexRegion.RegionType.HYPER, "foo[1..4)HX");
    check(ComplexRegion.RegionType.COMPLEX, "foo[1..4)CX");
    check(ComplexRegion.RegionType.OVERCOVERAGE, "foo[1..4)OC");
    check(ComplexRegion.RegionType.COMPLEX_NO_VARIANT, "foo[1..4)CXF");
    check(ComplexRegion.RegionType.NO_HYPOTHESES, "foo[1..4)NH");
    check(ComplexRegion.RegionType.TOO_MANY_HYPOTHESES, "foo[1..4)TMH");
  }

  protected void check(final RegionType type, final String exp) {
    final ComplexRegion cr = new ComplexRegion("foo", 1, 4, type);
    cr.integrity();
    assertEquals(type, cr.type());
    assertEquals(4, cr.getEnd());
    assertEquals(1, cr.getStart());
    assertEquals(exp, cr.toString());
  }

  public void testIn() {
    final RegionType type = ComplexRegion.RegionType.INTERESTING;
    final ComplexRegion cr = new ComplexRegion("foo", 1, 2, type);
    cr.integrity();
    assertEquals(type, cr.type());
    assertEquals(2, cr.getEnd());
    assertEquals(1, cr.getStart());
    assertEquals("foo[1..2)IN", cr.toString());
  }


}
