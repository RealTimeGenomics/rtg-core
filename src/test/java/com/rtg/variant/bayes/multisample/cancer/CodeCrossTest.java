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

package com.rtg.variant.bayes.multisample.cancer;

import com.rtg.variant.bayes.Code;

import junit.framework.TestCase;

/**
 */
public class CodeCrossTest extends TestCase {

  public void test() {

  }

  public void testValid() {
    final Code tc = new CodeCross(4);
    assertEquals(16, tc.size());
    assertFalse(tc.valid(-1));
    assertFalse(tc.valid(16));
    for (int i = 0; i < 16; ++i) {
      assertTrue(tc.valid(i));
    }
  }

  private void checkab(final Code tc, final int n, final int a, final int b) {
    assertEquals(a, tc.a(n));
    assertEquals(b, tc.b(n));
    assertEquals(b, tc.bc(n));
  }

  public void test4ab() {
    final Code tc = new CodeCross(3);
    assertEquals(9, tc.size());
    checkab(tc, 0, 0, 0);
    checkab(tc, 1, 0, 1);
    checkab(tc, 2, 0, 2);
    checkab(tc, 3, 1, 0);
    checkab(tc, 4, 1, 1);
    checkab(tc, 5, 1, 2);
    checkab(tc, 6, 2, 0);
    checkab(tc, 7, 2, 1);
    checkab(tc, 8, 2, 2);
  }

  public void testCode2() {
    final Code tc = new CodeCross(4);
    assertEquals(4, tc.rangeSize());
    for (int i = 0; i < tc.size(); ++i) {
      final int a = tc.a(i);
      assertTrue(0 <= a && a < 4);
      final int b = tc.bc(i);
      assertTrue(0 <= b && b < 4);
      assertEquals(i, tc.code(a, b));
    }
  }

  public void testCode() {
    final Code tc = new CodeCross(3);
    for (int i = 0; i < 3; ++i) {
      assertEquals(i, tc.code(i));
      assertTrue(tc.homozygous(tc.code(i, i)));
      for (int j = 0; j < 3; ++j) {
        if (i == j) {
          continue;
        }
        assertFalse(tc.homozygous(tc.code(i, j)));
      }

    }
  }
}
