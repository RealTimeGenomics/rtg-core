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
    for (int k = 0; k < delta.length; ++k) {
      assert a[k].length == b[k].length;
      delta[k] = new double[a[k].length];
      for (int j = 0; j < delta[k].length; ++j) {
        delta[k][j] = arith.poss2Prob(a[k][j]) - arith.poss2Prob(b[k][j]);
      }
    }
    return frobeniusNorm(delta);
  }
}
