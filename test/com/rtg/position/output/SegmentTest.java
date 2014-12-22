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
