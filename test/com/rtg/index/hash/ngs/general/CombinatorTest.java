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


import junit.framework.TestCase;

/**
 */
public class CombinatorTest extends TestCase {

  public void test() {
    check(4, 2);
    check(0, 0);
    check(1, 0);
    check(1, 1);
    check(5, 3);
  }

  private void check(final int n, final int i) {
    final int[] count = new int[1];
    final Combinator comb = new Combinator(n, i) {
      @Override
      public void permutation() {
        count[0]++;
        int c = 0;
        for (int x = 0; x < size(); ++x) {
          if (getBit(x)) {
            ++c;
          }
        }
        assertEquals(toBeSet(), c);
      }
    };
    comb.combine();
    assertEquals(Util.binomial(n, i), count[0]);
  }
}
