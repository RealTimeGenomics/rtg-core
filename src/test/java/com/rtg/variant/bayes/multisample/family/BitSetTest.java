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

package com.rtg.variant.bayes.multisample.family;

import junit.framework.TestCase;

/**
 */
public class BitSetTest extends TestCase {

  public void test() {
    final BitSet bs = BitSet.DNA_SET;
    bs.integrity();
    assertEquals(4, bs.length());
    assertEquals("A:C:G:T", bs.toString());
    assertEquals("empty", bs.toString(0));
    assertEquals("A", bs.toString(1 << 0));
    assertEquals("C", bs.toString(1 << 1));
    assertEquals("G", bs.toString(1 << 2));
    assertEquals("T", bs.toString(1 << 3));

    assertEquals(0, bs.toSet());
    assertEquals("A", bs.toString(bs.toSet(0)));
    assertEquals("C", bs.toString(bs.toSet(1)));
    assertEquals("G", bs.toString(bs.toSet(2)));
    assertEquals("T", bs.toString(bs.toSet(3)));

    assertTrue(bs.contains(1, 0));
    assertFalse(bs.contains(1, 1));
    assertFalse(bs.contains(1, 2));
    assertFalse(bs.contains(1, 3));

    assertTrue(bs.contains(15, 0));
    assertTrue(bs.contains(15, 1));
    assertTrue(bs.contains(15, 2));
    assertTrue(bs.contains(15, 3));

    assertEquals(12, bs.complement(3));
    assertEquals(15, bs.complement(0));
  }

  public void test1() {
    final BitSet bs = new BitSet("A");
    bs.integrity();
    assertEquals(1, bs.length());
    assertEquals("A", bs.toString());
    assertEquals("empty", bs.toString(0));
    assertEquals("A", bs.toString(1 << 0));

    assertEquals(0, bs.toSet());
    assertEquals("A", bs.toString(bs.toSet(0)));

    assertEquals(0, bs.complement(1));
    assertEquals(1, bs.complement(0));
  }

  //32 entries - special case for bit handling
  public void test32() {
    final String[] s = new String[32];
    for (int i = 0; i < 32; ++i) {
      s[i] = "" + i;
    }
    final BitSet bs = new BitSet(s);
    bs.integrity();
    assertEquals(32, bs.length());
    assertEquals("empty", bs.toString(0));
    for (int i = 0; i < 32; ++i) {
      assertEquals(s[i], bs.toString(1 << i));
    }

    assertEquals(s[0] + ":" + s[31], bs.toString(bs.toSet(0, 31)));
    assertEquals(0, bs.toSet());

    assertEquals(Integer.MAX_VALUE, bs.complement(Integer.MIN_VALUE));
    assertEquals(-1, bs.complement(0));

  }
}
