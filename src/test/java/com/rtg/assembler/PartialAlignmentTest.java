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

package com.rtg.assembler;

import junit.framework.TestCase;

/**
 */
public class PartialAlignmentTest extends TestCase {
  public void test() {
    final PartialAlignment pa = new PartialAlignment(2, 10, 20, 30, 40, 50);
    assertFalse(pa.equals(null));
    assertFalse(pa.equals("ASDF"));
    assertTrue(pa.equals(pa));
    PartialAlignment pa2 = new PartialAlignment(3, 10, 20, 30, 40, 50);
    assertFalse(pa.equals(pa2));
    pa2 = new PartialAlignment(3, 11, 20, 30, 40, 50);
    assertFalse(pa.equals(pa2));
    pa2 = new PartialAlignment(2, 10, 22, 30, 40, 50);
    assertFalse(pa.equals(pa2));
    pa2 = new PartialAlignment(2, 10, 20, 32, 40, 50);
    assertFalse(pa.equals(pa2));
    pa2 = new PartialAlignment(2, 10, 20, 30, 42, 50);
    assertFalse(pa.equals(pa2));
    pa2 = new PartialAlignment(2, 10, 20, 30, 40, 52);
    assertFalse(pa.equals(pa2));
    pa2 = new PartialAlignment(2, 10, 20, 30, 40, 50);
    assertTrue(pa.equals(pa2));
    assertEquals(pa.hashCode(), pa2.hashCode());
    assertEquals(-916452864, pa.hashCode());
  }
  public void testGetters() {
    final PartialAlignment pa = new PartialAlignment(2, 10, 20, 30, 40, 50);
    assertEquals(2, pa.getAlignmentScore());
    assertEquals(10, pa.getReadStart());
    assertEquals(20, pa.getReadEnd());
    assertEquals(30, pa.getContig());
    assertEquals(40, pa.getContigStart());
    assertEquals(50, pa.getContigEnd());
    assertEquals("contig=30 [40,50] read [10,20] score=2", pa.toString());
  }
}
