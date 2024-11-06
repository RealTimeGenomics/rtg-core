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
