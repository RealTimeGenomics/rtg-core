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

import com.rtg.metagenomics.matrix.Vector;

import junit.framework.TestCase;

/**
 */
public class LLineTest extends TestCase {

  //see specieslinetest.xlsx
  public void test() {
    final Frag f0 = SpeciesLineTest.getFrag(new int[]{0}, new int[]{3}, 1);
    final Frag f1 = SpeciesLineTest.getFrag(new int[]{1}, new int[]{4}, 2);
    final Frag f2 = SpeciesLineTest.getFrag(new int[]{2}, new int[]{5}, 1);
    final Frag f02 = SpeciesLineTest.getFrag(new int[] {0, 2}, new int[] {1, 2}, 1);
    final Frag f12 = SpeciesLineTest.getFrag(new int[] {0, 1}, new int[] {2, 3}, 1);
    final Frag[] frags = {f0, f1, f2, f02, f12};
    final long[] lengths = {10, 15, 20};
    final BlockInfo bi = new BlockInfo(42, null , frags, null, lengths, false);
    final Vector r = new Vector(new double[] {0.1, 0.2, 0.3});
    final Vector d = new Vector(new double[] {-0.1, +0.1, -0.2});

    final Line spl = new LLine(r, d, bi);
    assertEquals(1, spl.derivativeOrder());
    final double[] vs0 = spl.values(0.0);
    assertEquals(2, vs0.length);

    assertEquals(11.82461, vs0[0], 1e-5);
    assertEquals(-0.76429, vs0[1], 1e-5);
    final double[] vs1 = spl.values(1.0);
    assertEquals(2, vs1.length);
    assertEquals(11.18881, vs1[0], 1e-5);
    assertEquals(-0.51408, vs1[1], 1e-5);
  }

}
