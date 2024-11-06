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

package com.rtg.index;

import java.util.Arrays;

import com.rtg.AbstractTest;

/**
 * Test
 */
public class UnionRepeatFrequencyFilterMethodTest extends AbstractTest {

  public void test() {
    final IndexFilterMethod m = new UnionRepeatFrequencyFilterMethod(
      new FixedRepeatFrequencyFilterMethod(50),
      new BlacklistFilterMethod(Arrays.asList(0b1111L, 0b1010L, 0b1100L), 4, 1));
    m.initialize(null);
    // Fail blacklist, (freq OK)
    assertFalse(m.keepHash(0b1111L, 10));
    assertFalse(m.keepHash(0b1100L, 10));
    assertFalse(m.keepHash(0b1010L, 10));
    // Fail frequency, (blacklist OK)
    assertFalse(m.keepHash(0, 80));
    assertFalse(m.keepHash(0, 1000));
    assertFalse(m.keepHash(0, 10000));
    assertFalse(m.keepHash(0b1011L, 60));
    // Pass both
    assertTrue(m.keepHash(0b1110L, 10));
    assertTrue(m.keepHash(0b1101L, 10));
    assertTrue(m.keepHash(0b1011L, 10));
    assertTrue(m.keepHash(0, 5));
    assertTrue(m.keepHash(0, 25));
    assertTrue(m.keepHash(0, 50));
  }
}
