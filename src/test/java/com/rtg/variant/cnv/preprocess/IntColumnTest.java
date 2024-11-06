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

package com.rtg.variant.cnv.preprocess;

import junit.framework.TestCase;

/**
 * Tests the corresponding class.
 */
public class IntColumnTest extends TestCase {

  public void test() {
    final IntColumn col = new IntColumn("col");
    assertEquals("col", col.getName());
    col.add(42);
    assertEquals(1, col.size());
    assertEquals(42.0, col.get(0));
    assertEquals("42", col.toString(0));
    assertEquals(42.0, col.mean());
    assertEquals(42.0, col.median());
    col.remove(0);
    assertEquals(0, col.size());
    col.add(1);
    col.add(2);
    col.add(3);
    col.add(4);
    assertEquals(4, col.size());
    col.remove(1);  // 1,3,4
    assertEquals(3, col.size());
    assertEquals(1, (int) col.get(0));
    assertEquals(3, (int) col.get(1));
    col.add(0, 3);  // 3,1,3,4
    assertEquals(3, (int) col.get(0));
    assertEquals(1, (int) col.get(1));
    assertEquals(3, (int) col.get(2));
    assertEquals(4, (int) col.get(3));
    col.add(3, 5);  // 3,1,3,5,4
    assertEquals(5, col.size());
    assertEquals(4, (int) col.get(4));
  }
}
