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

import junit.framework.TestCase;



/**
 * Test SegmentCollection
 */
public class SegmentCollectionTest extends TestCase {

  /**
   * Check clear()
   */
  public void test1() {
    final SegmentCollection sc = new SegmentCollection(10);
    sc.globalIntegrity();
    assertEquals(0, sc.size());
    check(sc, 0, "0:0");
    check(sc, -1, "-1:0");
    assertEquals("SegmentCollection [0]" + StringUtils.LS, sc.toString());

    final Segment s1 = new Segment();
    s1.initialize(2, 0, 1);
    sc.add(s1);
    sc.globalIntegrity();
    assertEquals(1, sc.size());
    assertEquals(s1, sc.get(0));
    check(sc, 1, "1:1");
    check(sc, -1, "-1:1");
    assertEquals("SegmentCollection [1]" + StringUtils.LS + "[0] 2:0..1" + StringUtils.LS, sc.toString());

    sc.clear();
    sc.globalIntegrity();
    assertEquals(0, sc.size());
    check(sc, 0, "0:0");
    check(sc, -1, "-1:0");
    assertEquals("SegmentCollection [0]" + StringUtils.LS, sc.toString());
  }

  /**
   * Check removeNext()
   */
  public void test2() {
    final SegmentCollection sc = new SegmentCollection(10);
    sc.globalIntegrity();
    assertEquals(0, sc.size());
    assertEquals(null, sc.removeNext());
    check(sc, 0, "0:0");
    check(sc, -1, "-1:0");
    assertEquals("SegmentCollection [0]" + StringUtils.LS, sc.toString());

    final Segment s1 = new Segment();
    s1.initialize(2, 0, 1);
    sc.add(s1);
    sc.globalIntegrity();
    assertEquals(1, sc.size());
    assertEquals(s1, sc.get(0));
    check(sc, 1, "1:1");
    check(sc, -1, "-1:1");
    assertEquals("SegmentCollection [1]" + StringUtils.LS + "[0] 2:0..1" + StringUtils.LS, sc.toString());

    assertEquals(s1, sc.removeNext());
    sc.globalIntegrity();
    assertEquals(0, sc.size());
    assertEquals(null, sc.removeNext());
    check(sc, 0, "0:0");
    check(sc, -1, "-1:0");
    assertEquals("SegmentCollection [0]" + StringUtils.LS, sc.toString());
  }

  /**
   * Check flush()
   * @throws IOException
   */
  public void test3() throws IOException {
    final SegmentCollection sc = new SegmentCollection(10);
    sc.globalIntegrity();
    assertEquals(0, sc.size());
    assertEquals(null, sc.removeNext());
    check(sc, 0, "0:0");
    check(sc, -1, "-1:0");
    assertEquals("SegmentCollection [0]" + StringUtils.LS, sc.toString());

    final Segment s1 = new Segment();
    s1.initialize(2, 0, 1);
    sc.add(s1);
    sc.globalIntegrity();
    assertEquals(1, sc.size());
    assertEquals(s1, sc.get(0));
    check(sc, 1, "1:1");
    check(sc, -1, "-1:1");
    assertEquals("SegmentCollection [1]" + StringUtils.LS + "[0] 2:0..1" + StringUtils.LS, sc.toString());

    final Segment s2 = new Segment();
    s2.initialize(2, 0, 12);
    sc.add(s2);
    sc.globalIntegrity();
    assertEquals(2, sc.size());
    assertEquals(s1, sc.get(0));
    assertEquals(s2, sc.get(1));
    check(sc, 2, "2:2");
    check(sc, -1, "-1:2");
    assertEquals("SegmentCollection [2]" + StringUtils.LS + "[0] 2:0..1" + StringUtils.LS + "[1] 2:0..12" + StringUtils.LS, sc.toString());

    final SegmentCollection free = new SegmentCollection(2);
    free.globalIntegrity();
    final SegmentWriter out = new MockSegmentWriter();
    sc.flush(out, free, 42);
    sc.globalIntegrity();
    assertEquals(0, sc.size());
    free.globalIntegrity();
    assertEquals(2, free.size());
    assertEquals("2\t0\t2" + StringUtils.LS + "2\t0\t13" + StringUtils.LS, out.toString());
  }

  private static void check(final SegmentCollection sc, final int index, final String msg) {
    try {
      sc.get(index);
      fail();
    } catch (final IndexOutOfBoundsException e) {
      //expected
      assertEquals(msg, e.getMessage());
    }
  }
}
