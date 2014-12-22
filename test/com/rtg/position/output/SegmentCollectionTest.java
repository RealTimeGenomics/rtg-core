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
