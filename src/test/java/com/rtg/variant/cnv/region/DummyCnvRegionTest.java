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
