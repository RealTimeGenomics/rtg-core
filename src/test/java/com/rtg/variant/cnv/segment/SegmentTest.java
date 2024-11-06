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

package com.rtg.variant.cnv.segment;

import junit.framework.TestCase;

/**
 * Tests the corresponding class.
 */
public class SegmentTest extends TestCase {

  public void test() {
    final Segment left = new Segment("test", 42, 43, 100, 1, 1, 1);
    assertEquals(1, left.bins());
    assertEquals(100.0, left.sum());
    assertEquals(10000.0, left.sumSquares());
    assertEquals(100.0, left.mean());
    assertEquals(0.0, left.meanDistanceBetween());
    assertEquals(42, left.getStart());
    assertEquals(43, left.getEnd());
    final Segment right = new Segment("test", 43, 45, 50, 1, 1, 1);
    final Segment m = new Segment(left, right, Math.PI);
    assertEquals(left, m.left());
    assertEquals(right, m.right());
    assertEquals(Math.PI, m.deltaEnergy());
    assertEquals(2, m.bins());
    assertEquals(150.0, m.sum());
    assertEquals(12500.0, m.sumSquares());
    assertEquals(75.0, m.mean());
    assertEquals(1.0, m.meanDistanceBetween());
    assertEquals(42, m.getStart());
    assertEquals(45, m.getEnd());
    assertEquals("2", m.toString());
    assertEquals(1.0, m.distanceToPrevious());
    assertEquals(1, m.firstBinLength());
    assertEquals(2, m.lastBinLength());
  }

  public void testAbsorbLeft() {
    final Segment left = new Segment("test", 42, 43, 100, 1, 1, 1);
    final Segment rl = new Segment("test", 43, 44, 25, 1, 1, 1);
    final Segment rr = new Segment("test", 44, 45, 25, 1, 1, 1);
    final Segment right = new Segment(rl, rr, 1.0);
    final Segment seg = Segment.absorbLeft(left, right);
    assertEquals(2, seg.bins());
    assertEquals(42, seg.getStart());
    assertEquals(45, seg.getEnd());
    assertEquals("test", seg.getSequenceName());
    assertEquals(50.0, seg.sum());
    assertEquals(1, seg.left().bins());
    assertEquals(1, seg.right().bins());
  }

  public void testAbsorbRight() {
    final Segment ll = new Segment("test", 42, 43, 50, 1, 1, 1);
    final Segment lr = new Segment("test", 43, 44, 25, 1, 1, 1);
    final Segment left = new Segment(ll, lr, Math.PI);
    final Segment right = new Segment("test", 44, 45, 100, 1, 1, 1);
    final Segment seg = Segment.absorbRight(left, right);
    assertEquals(2, seg.bins());
    assertEquals(42, seg.getStart());
    assertEquals(45, seg.getEnd());
    assertEquals("test", seg.getSequenceName());
    assertEquals(75.0, seg.sum());
    assertEquals(Math.PI, seg.deltaEnergy());
    assertEquals(1, seg.left().bins());
    assertEquals(1, seg.right().bins());
  }
}
