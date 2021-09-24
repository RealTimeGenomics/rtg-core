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

import junit.framework.TestCase;

/**
 */
public class BinomialSpecialTest extends TestCase {

  public void test() {
    check(1, 1);
    check(2, 1);
    check(4, 2);
    check(4, 4);
    check(5, 5);
    check(6, 5);
    check(10, 5);
    check(100, 5);
    check(200, 5);
    check(BinomialSpecial.LENGTH - 2, 5);
    check(BinomialSpecial.LENGTH - 1, 5);
    check(BinomialSpecial.LENGTH, 5);
    check(BinomialSpecial.LENGTH + 1, 5);
  }

  private void check(int n, int a) {
    final double exp = MathUtils.logBinomial(n, a);
    final double actual = BinomialSpecial.logBinomial(n, a);
    //System.err.println(n + ":" + a + " " + Math.exp(exp) + ":" + Math.exp(actual));
    assertEquals(exp, actual, 1e-7);
  }
}
