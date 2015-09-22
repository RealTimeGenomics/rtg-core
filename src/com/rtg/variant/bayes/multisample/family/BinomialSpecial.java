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
    for (int i = 1; i < LENGTH; i++) {
      PASCAL[0][i] = 1;
      for (int j = 1; j < A_LIMIT; j++) {
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
