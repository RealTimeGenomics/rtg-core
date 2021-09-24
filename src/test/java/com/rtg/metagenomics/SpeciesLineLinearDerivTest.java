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
