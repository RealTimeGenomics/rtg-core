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

package com.rtg.segregation;

import junit.framework.TestCase;

/**
 */
public class CrossOverTest extends TestCase {

  public void test0() {
    final PatternArray pa = new PatternArray("0", "1");
    final CrossOver co = new CrossOver(3, false, 42, 2, 1, pa);
    co.integrity();
    assertTrue(pa == co.pattern());
    assertEquals("mo 42 2 1", co.toString());

    final int[] bedSearch0 = co.searchLevel(3, 2);
    assertEquals(3, bedSearch0[0]);
    assertEquals(1, bedSearch0[1]);

    final int[] bedSearch1 = co.searchLevel(3, 1);
    assertEquals(3, bedSearch1[0]);
    assertEquals(2, bedSearch1[1]);
  }

  public void test1() {
    final PatternArray pa = new PatternArray("0", "1");
    final CrossOver co = new CrossOver(3, true, 0, 1, 2, pa);
    co.integrity();
    assertTrue(pa == co.pattern());
    assertEquals("fa 0 1 2", co.toString());

    final int[] bedSearch0 = co.searchLevel(1, 3);
    assertEquals(2, bedSearch0[0]);
    assertEquals(3, bedSearch0[1]);

    final int[] bedSearch1 = co.searchLevel(2, 3);
    assertEquals(1, bedSearch1[0]);
    assertEquals(3, bedSearch1[1]);
  }
}
