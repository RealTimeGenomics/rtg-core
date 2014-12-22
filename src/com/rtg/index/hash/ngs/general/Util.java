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
package com.rtg.index.hash.ngs.general;

/**
 */
public final class Util {

  private Util() { }

  private static final long[][] BINOMIAL = new long[65][];
  static {
    for (int k = 1; k < BINOMIAL.length; k++) {
      final long[] row = new long[k];
      for (int j = 0; j < row.length; j++) {
        row[j] = binomial(k - 1, j) + binomial(k - 1, j + 1);
      }
      BINOMIAL[k] = row;
    }
  }

  /**
   * Compute binomial coefficient for k out of n.
   * @param n total number.
   * @param k number selected from total.
   * @return the number of ways of selecting k from n.
   */
  public static long binomial(final int n, final int k) {
    if (k > n || k < 0) {
      return 0;
    }
    if (k == 0) {
      return 1;
    }
    return BINOMIAL[n][k - 1];
  }
}
