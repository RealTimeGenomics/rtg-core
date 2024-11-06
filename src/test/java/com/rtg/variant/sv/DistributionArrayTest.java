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
