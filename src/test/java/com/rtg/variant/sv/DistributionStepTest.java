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
