/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
