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
public class SpeciesLineLinearDerivTest extends TestCase {

  //see specieslinetest.xlsx
  public void test() {
    final Frag[] frags = SpeciesLineLinearTest.buildFrags();
    final long[] lengths = {10, 15, 20};
    final BlockInfo bi = new BlockInfo(42, null , frags, null, lengths, false);
    final Vector r = new Vector(new double[] {0.1, 0.2, 0.3});
    final Vector d = new Vector(new double[] {-0.1, +0.1, -0.2});
    final Line spl = new SpeciesLineLinearDeriv(r, d, bi);
    assertEquals(1, spl.derivativeOrder());
    final double val = -2.24405;
    assertEquals(val, spl.value(0.0), 1e-5);
    final double[] values = spl.values(0.0);
    assertEquals(2, values.length);
    assertEquals(val, values[0], 1e-5);
    assertEquals(2.47027, values[1], 1e-5);


    final double[] val1 = spl.values(0.5);
    assertEquals(2, val1.length);
    assertEquals(-0.30654, val1[0], 1e-5);
    assertEquals(6.56841, val1[1], 1e-5);
  }

}
