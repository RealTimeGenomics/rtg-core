/*
 * Copyright (c) 2014. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
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
