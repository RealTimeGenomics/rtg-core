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
public class SortedMultiSetTest extends TestCase {
  public void test() {
    final SortedMultiSet<Integer> set = new SortedMultiSet<>();
    assertEquals("[]", set.toString());
    assertEquals("{}", set.countMap().toString());
    assertEquals(0, set.get(0));

    set.add(0);
    assertEquals("[0->1]", set.toString());
    assertEquals("{0=1}", set.countMap().toString());
    assertEquals(1, set.get(0));

    set.add(0);
    assertEquals("[0->2]", set.toString());
    assertEquals("{0=2}", set.countMap().toString());
    assertEquals(2, set.get(0));

    for (int i = 1; i < 10; ++i) {
      set.add(i);
    }
    final String exp10 = "[0->2, 1->1, 2->1, 3->1, 4->1, 5->1, 6->1, 7->1, 8->1, 9->1";
    assertEquals(exp10 + "]", set.toString());
    assertEquals(2, set.get(0));
    for (int i = 1; i < 10; ++i) {
      assertEquals(1, set.get(i));
    }

    set.add(12);
    assertEquals(exp10 + StringUtils.LS + ", 12->1"  + StringUtils.LS + "]", set.toString());
    assertEquals(2, set.get(0));
    for (int i = 1; i < 10; ++i) {
      assertEquals(1, set.get(i));
    }
    assertEquals(0, set.get(11));
    assertEquals(1, set.get(12));

  }
}
