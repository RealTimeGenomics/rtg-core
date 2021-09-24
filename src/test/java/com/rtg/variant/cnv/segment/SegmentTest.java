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
