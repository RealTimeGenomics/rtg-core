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
package com.rtg.position.output;

import com.rtg.util.TestUtils;

import junit.framework.TestCase;

/**
 * Test <code>Segment</code>.
 */
public class SegmentTest extends TestCase {

  public void test1() {
    final Segment s = new Segment();
    s.integrity();
    assertTrue(s.isEmpty());
    assertEquals("empty", s.toString());

    s.initialize(2, 0, 1);
    s.integrity();
    assertFalse(s.isEmpty());
    assertEquals(2, s.seqId());
    assertEquals(0, s.start());
    assertEquals(1, s.end());
    assertEquals("2:0..1", s.toString());

    s.extend(3);
    s.integrity();
    assertFalse(s.isEmpty());
    assertEquals(2, s.seqId());
    assertEquals(0, s.start());
    assertEquals(3, s.end());
    assertEquals("2:0..3", s.toString());

    s.clear();
    s.integrity();
    assertTrue(s.isEmpty());
    assertEquals("empty", s.toString());

    //check case when start != 0
    s.initialize(2, 7, 9);
    s.integrity();
    assertFalse(s.isEmpty());
    assertEquals(2, s.seqId());
    assertEquals(7, s.start());
    assertEquals(9, s.end());
    assertEquals("2:7..9", s.toString());
  }

  /**
   * There was a problem in the wild with length 1 sequences.
   */
  public void testLength1() {
    final Segment s = new Segment();
    s.integrity();
    assertTrue(s.isEmpty());
    assertEquals("empty", s.toString());

    s.initialize(2, 0, 0);
    s.integrity();
    assertFalse(s.isEmpty());
    assertEquals(2, s.seqId());
    assertEquals(0, s.start());
    assertEquals(0, s.end());
    assertEquals("2:0..0", s.toString());

    s.clear();
    s.initialize(2, 1, 1);
    s.integrity();
    assertFalse(s.isEmpty());
    assertEquals(2, s.seqId());
    assertEquals(1, s.start());
    assertEquals(1, s.end());
    assertEquals("2:1..1", s.toString());
  }

  public void testOrder() {
    final Segment s00 = init(0, 0, 1);
    final Segment s01 = init(0, 0, 1);
    final Segment s02 = init(0, 0, 2);
    final Segment s10 = init(1, 0, 1);
    final Segment s1i = init(1, Integer.MAX_VALUE - 1, Integer.MAX_VALUE);
    final Segment si0 = init(Integer.MAX_VALUE, 1, 2);
    TestUtils.equalsHashTest(new Segment[][] {{s00, s01}, {s02}, {s10}, {s1i}, {si0}});
    TestUtils.testOrder(new Segment[] {s00, s02, s10, s1i, si0}, true);
    assertEquals(0, s00.compareTo(s01));
   }

  private Segment init(final int a, final int b, final int c) {
    final Segment s = new Segment();
    s.integrity();
    s.initialize(a, b, c);
    s.integrity();
    return s;
  }
}
