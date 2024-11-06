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

package com.rtg.variant.bayes.multisample.family;

import com.rtg.util.MathUtils;

/**
 * Quickly compute ln (binomial(n, m)) where 1 &le; m &le; 5.
 */
public final class BinomialSpecial {

  private static final int A_LIMIT = 6;

  private BinomialSpecial() { }

  /** Determined by running main and finding first position that overflows an int. */
  static final int LENGTH = 16175;

  // Uses a table up to n = A_LIMIT which allows the table to be computed using Pascals triangle.
  // Otherwise it uses logFactorial which is slow but robust.
  private static final long[][] PASCAL = new long[A_LIMIT][LENGTH];
  private static final double[][] PASCAL_LN = new double[A_LIMIT][LENGTH];
  static {
    PASCAL[0][0] = 1;
    for (int i = 1; i < LENGTH; ++i) {
      PASCAL[0][i] = 1;
      for (int j = 1; j < A_LIMIT; ++j) {
        PASCAL[j][i] = PASCAL[j - 1][i - 1] + PASCAL[j][i - 1];
        PASCAL_LN[j][i] = Math.log(PASCAL[j][i]);
      }
    }
  }

  /**
   * Compute <code>ln(binomial(n, a))</code> quickly.
   * @param n total count.
   * @param a sub-count.
   * @return <code>ln(binomial(n, a))</code>.
   */
  public static double logBinomial(final int n, final int a) {
    assert 1 <= a && a < A_LIMIT && a <= n;
    if (n < LENGTH) {
      return PASCAL_LN[a][n];
    }
    //fall back to the slow way
    return MathUtils.logFactorial(n) - MathUtils.logFactorial(a) - MathUtils.logFactorial(n - a);
  }
}
