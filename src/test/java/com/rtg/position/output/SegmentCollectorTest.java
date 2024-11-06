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

import java.io.IOException;

import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;

import junit.framework.TestCase;



/**
 */
public class SegmentCollectorTest extends TestCase {

  /**
   * Tests most cases in the loop inside endPosition.
   * Ensures the end case when the values in the collector are exhausted first is covered.
   * @throws IOException
   */
  public void test1() throws IOException {
    final MockSegmentWriter wr = new MockSegmentWriter();
    final SegmentCollector sc = new SegmentCollector(3, 3, 1, wr);
    sc.globalIntegrity();
    assertEquals(0, sc.size());

    sc.add(0, 5);
    sc.globalIntegrity();
    assertEquals(1, sc.size());

    sc.add(0, 3);
    sc.globalIntegrity();
    assertEquals(2, sc.size());

    sc.add(0, 25);
    sc.globalIntegrity();
    assertEquals(3, sc.size());

    final SegmentCollection s0 = new SegmentCollection(4);
    add(s0, 0, 0, 2);
    add(s0, 0, 20, 22);
    add(s0, 1, 0, 2);
    add(s0, 1, 0, 5);

    final SegmentCollection free = new SegmentCollection(2);
    final Segment sf = new Segment();
    free.add(sf);

    final SegmentCollection newS = new SegmentCollection(3);
    sc.endPosition(s0, newS, free, 42);
    assertEquals(0, s0.size());
    s0.globalIntegrity();
    assertEquals(3, newS.size());
    final String xNewS = ""
    + "SegmentCollection [3]" + StringUtils.LS
    + "[0] 0:0..3" + StringUtils.LS
    + "[1] 0:3..5" + StringUtils.LS
    + "[2] 0:23..25" + StringUtils.LS
    ;
    assertEquals(xNewS, newS.toString());
    newS.globalIntegrity();

    assertEquals(2, free.size());
    free.globalIntegrity();

    final String xWr = ""
      + "0\t20\t3" + StringUtils.LS
      + "1\t0\t3" + StringUtils.LS
      + "1\t0\t6" + StringUtils.LS
      ;
    assertEquals(xWr, wr.toString());
  }

  /**
   * Tests most cases in the loop inside endPosition.
   * Ensures the end case when the oldS is exhausted first is covered.
   * @throws IOException
   */
  public void test2() throws IOException {
    final MockSegmentWriter wr = new MockSegmentWriter();
    final SegmentCollector sc = new SegmentCollector(4, 3, 1, wr);
    sc.globalIntegrity();
    assertEquals(0, sc.size());

    sc.add(0, 5);
    sc.globalIntegrity();
    assertEquals(1, sc.size());

    sc.add(1, 30);
    sc.globalIntegrity();
    assertEquals(2, sc.size());

    sc.add(0, 25);
    sc.globalIntegrity();
    assertEquals(3, sc.size());

    sc.add(0, 3);
    sc.globalIntegrity();
    assertEquals(4, sc.size());

    final SegmentCollection s0 = new SegmentCollection(2);
    add(s0, 0, 0, 2);
    add(s0, 0, 20, 22);

    final SegmentCollection free = new SegmentCollection(3);
    final Segment sf0 = new Segment();
    free.add(sf0);
    final Segment sf1 = new Segment();
    free.add(sf1);

    final SegmentCollection newS = new SegmentCollection(4);
    sc.endPosition(s0, newS, free, 42);
    assertEquals(4, newS.size());

    assertEquals(0, s0.size());
    s0.globalIntegrity();
    assertEquals(0, sc.size());

    final String xNewS = ""
    + "SegmentCollection [4]" + StringUtils.LS
    + "[0] 0:0..3" + StringUtils.LS
    + "[1] 0:3..5" + StringUtils.LS
    + "[2] 0:23..25" + StringUtils.LS
    + "[3] 1:28..30" + StringUtils.LS
    ;
    assertEquals(xNewS, newS.toString());
    newS.globalIntegrity();

    assertEquals(0, free.size());
    free.globalIntegrity();

    final String xWr = ""
      + "0\t20\t3" + StringUtils.LS
      ;
    assertEquals(xWr, wr.toString());
  }

  /**
   * Error message for exhausting free list.
   * @throws IOException
   */
  public void test3() throws IOException {
    final MockSegmentWriter wr = new MockSegmentWriter();
    final SegmentCollector sc = new SegmentCollector(4, 3, 1, wr);
    sc.globalIntegrity();
    assertEquals(0, sc.size());

    sc.add(0, 5);
    sc.globalIntegrity();
    assertEquals(1, sc.size());

    final SegmentCollection s0 = new SegmentCollection(2);

    final SegmentCollection free = new SegmentCollection(3);

    final SegmentCollection newS = new SegmentCollection(4);
    try {
      sc.endPosition(s0, newS, free, 42);
      fail();
    } catch (final RuntimeException e) {
      assertEquals("free list empty", e.getMessage());
    }
  }

  /**
   * Tests most cases in the loop inside endPosition.
   * Deal with case where seqId agrees but positions differ.
   * @throws IOException
   */
  public void test4() throws IOException {
    final MockSegmentWriter wr = new MockSegmentWriter();
    final SegmentCollector sc = new SegmentCollector(4, 3, 1, wr);
    sc.globalIntegrity();
    assertEquals(0, sc.size());

    sc.add(0, 5);
    sc.globalIntegrity();
    assertEquals(1, sc.size());

    sc.add(0, 30);
    sc.globalIntegrity();
    assertEquals(2, sc.size());

    sc.add(0, 25);
    sc.globalIntegrity();
    assertEquals(3, sc.size());

    sc.add(0, 3);
    sc.globalIntegrity();
    assertEquals(4, sc.size());

    final SegmentCollection s0 = new SegmentCollection(3);
    add(s0, 0, 0, 2);
    add(s0, 0, 20, 22);
    add(s0, 0, 40, 42);

    final SegmentCollection free = new SegmentCollection(3);
    final Segment sf0 = new Segment();
    free.add(sf0);
    final Segment sf1 = new Segment();
    free.add(sf1);

    final SegmentCollection newS = new SegmentCollection(4);
    sc.endPosition(s0, newS, free, 42);
    assertEquals(4, newS.size());

    assertEquals(0, s0.size());
    s0.globalIntegrity();
    assertEquals(0, sc.size());

    final String xNewS = ""
    + "SegmentCollection [4]" + StringUtils.LS
    + "[0] 0:0..3" + StringUtils.LS
    + "[1] 0:3..5" + StringUtils.LS
    + "[2] 0:23..25" + StringUtils.LS
    + "[3] 0:28..30" + StringUtils.LS
    ;
    assertEquals(xNewS, newS.toString());
    newS.globalIntegrity();

    assertEquals(1, free.size());
    free.globalIntegrity();

    final String xWr = ""
      + "0\t20\t3" + StringUtils.LS
      + "0\t40\t3" + StringUtils.LS
      ;
    assertEquals(xWr, wr.toString());
  }


  private void add(final SegmentCollection s0, final int seq, final int st, final int end) {
    final Segment g0 = new Segment();
    g0.initialize(seq, st, end);
    s0.add(g0);
  }

  public void testCode() {
    checkCode(0, 0);
    checkCode(0, 1);
    checkCode(1, 0);
    checkCode(Integer.MAX_VALUE, Integer.MAX_VALUE);
  }

  private void checkCode(final int seqId, final int posn) {
    final long v = SegmentCollector.pack(seqId, posn);
    assertEquals(seqId, SegmentCollector.seqId(v));
    assertEquals(posn, SegmentCollector.position(v));
  }

  public void testCodeOrder() {
    final long v0 = SegmentCollector.pack(0, 0);
    final long v1 = SegmentCollector.pack(0, 1);
    final long v2 = SegmentCollector.pack(0, 2);
    final long v3 = SegmentCollector.pack(0, Integer.MAX_VALUE);
    final long v4 = SegmentCollector.pack(1, 0);
    final long v5 = SegmentCollector.pack(1, Integer.MAX_VALUE);
    final long v6 = SegmentCollector.pack(Integer.MAX_VALUE, 0);
    final long v7 = SegmentCollector.pack(Integer.MAX_VALUE, Integer.MAX_VALUE);
    final Long[] sorted = {v0, v1, v2, v3, v4, v5, v6, v7};
    TestUtils.testOrder(sorted, false);
  }
}
