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

import java.util.Arrays;
import java.util.List;

import junit.framework.TestCase;

/**
 * Tests the corresponding class.
 */
public class SegmentChainTest extends TestCase {

  public void testNiceClean() {
    final SegmentChain sc = new SegmentChain(new EnergySegmentScorer(0, 0));
    sc.add(new Segment("test", 0, 1, 100, 0.0, 1, 1));
    sc.add(new Segment("test", 1, 2, 100, 0.0, 1, 1));
    sc.add(new Segment("test", 2, 3, 100, 0.0, 1, 1));
    sc.add(new Segment("test", 3, 4, 1, 0.0, 1, 1));
    sc.add(new Segment("test", 4, 5, 1, 0.0, 1, 1));
    sc.add(new Segment("test", 6, 7, 1, 0.0, 1, 1));
    sc.collapse();
    final List<Segment> segs = Arrays.asList(sc.get(0).left(), sc.get(0).right());
    assertEquals(2, segs.size());
    final Segment a = segs.get(0);
    assertEquals(segs.toString(), 3, a.bins());
    assertEquals(100.0, a.mean(), 1e-8);
    final Segment b = segs.get(1);
    assertEquals(3, b.bins());
    assertEquals(1.0, b.mean(), 1e-8);
  }

  public void testMiddleEnergy() {
    final SegmentChain sc = new SegmentChain(new EnergySegmentScorer(0, 0));
    for (int k = 0; k < 10; ++k) {
      sc.add(new Segment("test", k, k + 1, 1, 0.0, 1, 1));
    }
    sc.add(new Segment("test", 11, 12, 4, 0.0, 1, 1));
    sc.add(new Segment("test", 12, 13, 3, 0.0, 1, 1));
    for (int k = 0; k < 10; ++k) {
      sc.add(new Segment("test", 14 + k, 15 + k, 1, 0.0, 1, 1));
    }
    sc.collapse();
    final List<Segment> segs = Arrays.asList(sc.get(0).left().left(), sc.get(0).left().right(), sc.get(0).right());
    assertEquals(3, segs.size());
    final Segment a = segs.get(0);
    assertEquals(segs.toString(), 10, a.bins());
    assertEquals(1.0, a.mean(), 1e-8);
    final Segment b = segs.get(1);
    assertEquals(2, b.bins());
    assertEquals((4 + 3) / 2.0, b.mean(), 1e-8);
    final Segment c = segs.get(2);
    assertEquals(10, c.bins());
    assertEquals(1.0, c.mean(), 1e-8);
  }
}
