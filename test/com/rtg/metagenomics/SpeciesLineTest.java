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
import com.rtg.metagenomics.matrix.VectorSimple;

import junit.framework.TestCase;

/**
 */
public class SpeciesLineTest extends TestCase {

  static Frag getFrag(int[] genomes, int[] counts, int mult) {
    assert counts.length == genomes.length;
    final int n = genomes.length;
    int cnt = 0;
    for (int i = 0; i < n; i++) {
      cnt += counts[i];
    }
    return new Frag(genomes, counts, cnt, n, mult);
  }


  //see specieslinetest.xlsx
  public void test() {
    final Frag f0 = getFrag(new int[] {0}, new int[] {3}, 1);
    final Frag f1 = getFrag(new int[] {1}, new int[] {4}, 2);
    final Frag f2 = getFrag(new int[] {2}, new int[] {5}, 1);
    final Frag f02 = getFrag(new int[] {0, 2}, new int[] {1, 2}, 1);
    final Frag f12 = getFrag(new int[] {0, 1}, new int[] {2, 3}, 1);
    final Frag[] frags = {f0, f1, f2, f02, f12};
    final long[] lengths = {10, 15, 20};
    final BlockInfo bi = new BlockInfo(42, null , frags, null, lengths, false);
    final Vector r = new VectorSimple(new double[] {0.1, 0.2, 0.3});
    final Vector d = new VectorSimple(new double[] {-0.1, +0.1, -0.2});
    final Line spl = new SpeciesLine(r, d, bi);
    assertEquals(0, spl.derivativeOrder());
    final double val = -0.76428;
    assertEquals(val, spl.value(0.0), 1e-5);
    assertEquals(val, spl.values(0.0)[0], 1e-5);

    assertEquals(-0.51408, spl.value(1.0), 1e-5);
  }
}
