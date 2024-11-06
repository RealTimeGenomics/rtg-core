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

package com.rtg.variant.bayes;

import junit.framework.TestCase;

/**
 */
public class CodeHaploidTest extends TestCase {

  public void testHaploid() {
    final Code tc = new CodeHaploid(4);
    assertEquals(4, tc.rangeSize());
    assertEquals(4, tc.size());
    for (int i = 0; i < 4; ++i) {
      assertTrue(tc.homozygous(i));
      assertEquals(tc.a(i), tc.bc(i));
    }
  }

  public void testValid() {
    final Code tc = new CodeHaploid(4);
    assertFalse(tc.valid(-1));
    assertFalse(tc.valid(4));
    for (int i = 0; i < 4; ++i) {
      assertTrue(tc.valid(i));
      assertEquals(i, tc.a(i));
      assertEquals(i, tc.code(i));
      assertEquals(i, tc.code(i, i));
    }

    try {
      tc.b(0);
      fail();
    } catch (final UnsupportedOperationException e) {
      //expected
    }
    try {
      tc.code(-1);
      fail();
    } catch (final IllegalArgumentException e) {
      assertEquals("a=-1", e.getMessage());
    }
    try {
      tc.code(1, 2);
      fail();
    } catch (final UnsupportedOperationException e) {
      //expected
    }
  }
}
