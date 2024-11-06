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
package com.rtg.util;

import junit.framework.TestCase;

/**
 */
public class SizeSplitTest extends TestCase {

  public void test() {
    final SizeSplit ss = new SizeSplit(16, 3);
    assertEquals("SizeSplit: n=16 d=3 s=5 r=1", ss.toString());
    assertEquals(0, ss.start(0));
    assertEquals(6, ss.start(1));
    assertEquals(11, ss.start(2));
    assertEquals(16, ss.start(3));
  }

  public void test1() {
    check(16, 3);
    check(0, 100);
    check(100, 100);
    check(100, 99);
    check(100, 101);
  }

  void check(final int n, final int d) {
    final SizeSplit ss = new SizeSplit(n, d);
    final int s = n / d;
    int last = 0;
    for (int i = 1; i <= d; ++i) {
      final int st = ss.start(i);
      final int diff = st - ss.start(i - 1);
      assertTrue(last <= st);
      assertTrue(diff == s || diff == s + 1);
      last = st;
    }
    assertEquals(0, ss.start(0));
    assertEquals(n, ss.start(d));
  }
}
