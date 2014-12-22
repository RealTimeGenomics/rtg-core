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
package com.rtg.metagenomics.metasnp;

import com.rtg.variant.util.arithmetic.PossibilityArithmetic;

/**
 * Frobenius norm for a matrix.
 */
final class Frobenius {

  private Frobenius() { }

  static double frobeniusNorm(final double[][] a) {
    double s = 0;
    for (final double[] v : a) {
      for (final double e : v) {
        s += e * e;
      }
    }
    return Math.sqrt(s);
  }

  static double frobeniusDistance(final PossibilityArithmetic arith, final double[][] a, final double[][] b) {
    // Not expecting values to be way off from each other, so just go back
    // to ordinary arithmetic
    assert a.length == b.length;
    final double[][] delta = new double[a.length][];
    for (int k = 0; k < delta.length; k++) {
      assert a[k].length == b[k].length;
      delta[k] = new double[a[k].length];
      for (int j = 0; j < delta[k].length; j++) {
        delta[k][j] = arith.poss2Prob(a[k][j]) - arith.poss2Prob(b[k][j]);
      }
    }
    return frobeniusNorm(delta);
  }
}
