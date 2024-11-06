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
package com.rtg.metagenomics;

import java.io.IOException;
import java.util.ArrayList;

import com.rtg.launcher.AbstractNanoTest;
import com.rtg.util.integrity.Exam;

/**
 */
public class BlockMappingTest extends AbstractNanoTest {

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
