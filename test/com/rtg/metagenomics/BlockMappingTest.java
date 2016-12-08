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
package com.rtg.metagenomics;

import java.io.IOException;
import java.util.ArrayList;

import com.rtg.util.integrity.Exam;
import com.rtg.util.test.NanoRegression;

import junit.framework.TestCase;

/**
 */
public class BlockMappingTest extends TestCase {

  protected NanoRegression mNano;

  @Override
  public void setUp() {
    mNano = new NanoRegression(this.getClass());
  }

  @Override
  public void tearDown() throws IOException {
    try {
      mNano.finish();
    } finally {
      mNano = null;
    }
  }


  static Frag frag(final int...f) {
    final ArrayList<Integer> al = new ArrayList<>();
    for (final int fi : f) {
      al.add(fi);
    }
    return new Frag(al);
  }

  private void checkA(final int[] a, final int... b) {
    assertEquals(b.length, a.length);
    for (int i = 0; i < a.length; ++i) {
      assertEquals("" + i, a[i], b[i]);
    }
    for (int i = 1; i < a.length; ++i) {
      assertTrue("" + i, a[i] <= i);
    }
  }

  private void check(final int[] a, final int... b) {
    assertEquals(b.length, a.length);
    for (int i = 0; i < a.length; ++i) {
      assertEquals("" + i, a[i], b[i]);
    }
  }
  public void testConstructAComplicated() {
    final Frag[] frags = {
            frag(0, 1, 2),
            frag(3, 4),
            frag(5, 6),
            frag(3, 6),
            frag(1, 5),
            frag(4, 5)
    };
    final BlockMapping bm = new BlockMapping(frags, 7);
    bm.globalIntegrity();
    checkA(bm.getA(), 0, 0, 0, 0, 0, 0, 0);
  }

  public void testConstructAComplicated2() throws IOException {
    final Frag[] frags = {
            frag(0, 6),
            frag(5, 4),
            frag(4, 3),
            frag(3, 2),
    };
    final BlockMapping bm = new BlockMapping(frags, 7);
    bm.globalIntegrity();
    checkA(bm.getA(), 0, -1, 2, 2, 2, 2, 0);
    final int[][] b = bm.getBlocks();

    assertEquals(2, b.length);
    check(b[0], 0, 6);
    check(b[1], 2, 3, 4, 5);
    check(bm.getC(), 0, 0, 0, 1, 2, 3, 1);
    final String stats = bm.statistics();
    mNano.check("constructacomplicated2", stats);
  }

  public void testEmpty() throws IOException {
    final BlockMapping bm = new BlockMapping(new Frag[] {}, 2);
    mNano.check("empty", bm.statistics());
  }

  public void testFragStats() throws IOException {
    final BlockMapping bm = new BlockMapping(new Frag[] {frag(0, 1, 2), frag(3, 4), frag(5, 6)}, 7);
    mNano.check("fragstats", bm.statistics());
  }

  public void testConstructASimple() {
    final Frag fr = frag(1, 1, 42);
    final BlockMapping blockMapping = new BlockMapping(new Frag[]{fr}, 43);
    blockMapping.globalIntegrity();
    final int[] a = blockMapping.getA();

    Exam.assertEquals(-1, a[0]);
    Exam.assertEquals(1, a[1]);
    Exam.assertEquals(1, a[42]);
    //System.err.println(Arrays.toString(a));
  }
}
