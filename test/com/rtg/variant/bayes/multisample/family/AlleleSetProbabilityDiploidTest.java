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

import com.rtg.variant.bayes.Code;
import com.rtg.variant.bayes.CodeDiploid;

import junit.framework.TestCase;

/**
 */
public class AlleleSetProbabilityDiploidTest extends TestCase {

  public void testLookup() {
    final Code code = new CodeDiploid(5);
    check(1.0 * 1.0 * 1.0 / 1.0, code, 0, code.code(0, 0), code.code(0, 0));

    check(2.0 * 1.0 * 2.0 / 30.0, code, 0, code.code(0, 0), code.code(0, 1));
    check(1.0 * 1.0 * 1.0 / 30.0, code, 1, code.code(0, 0), code.code(0, 0));

    check(3.0 * 1.0 * 2.0 / 150.0, code, 0, code.code(0, 0), code.code(1, 2));
    check(1.0 * 1.0 * 1.0 / 150.0, code, 2, code.code(0, 0), code.code(1, 1));

    check(1.0 * 2.0 * 2.0 / 240.0, code, 3, code.code(0, 1), code.code(0, 2));
    check(4.0 * 2.0 * 2.0 / 240.0, code, 3, code.code(0, 1), code.code(2, 3));

    check(1.0 * 2.0 * 2.0 / 120.0, code, 4, code.code(0, 1), code.code(2, 3));
  }

  private void check(final double exp, final Code code, final int ref, final int f, final int m) {
    assertEquals(Math.log(exp), AlleleSetProbabilityDiploid.SINGLETON.probabilityLn(code, ref, f, m), 0.0);
    assertEquals(Math.log(exp), AlleleSetProbabilityDiploid.SINGLETON.probabilityLn(code, ref, m, f), 0.0);
  }

  //test symmetries and the two ways of accessing the table
  public void testLookupCode() {
    final Code code = new CodeDiploid(4);
    for (int i = 0; i < code.size(); ++i)  {
      for (int j = 0; j < code.size(); ++j)  {
        //System.err.println("i=" + i + " j=" + j);
        final byte ref = 1;
        final double exp = AlleleSetProbabilityDiploid.SINGLETON.probabilityLn(code, ref, i, j);
        assertEquals(exp, AlleleSetProbabilityDiploid.SINGLETON.probabilityLn(code, ref, j, i));
      }
    }
  }

  public void testG() {
    final double t = AlleleSetProbabilityDiploid.gSum();
    assertEquals(Math.pow(4, 5), t);
  }

  static int lCount(final double[][][][] lookup) {
    int cnt = 0;
    final int il0 = lookup[0].length;
    for (double[][][] aLookup : lookup) {
      final int il = aLookup.length;
      assertEquals(il0, il);
      final int kl0 = aLookup[0].length;
      for (final double[][] anALookup : aLookup) {
        final int kl = anALookup.length;
        assertEquals(kl0, kl);
        final int ll0 = anALookup[0].length;
        for (final double[] anAnALookup : anALookup) {
          final int ll = anAnALookup.length;
          assertEquals(ll0, ll);
          cnt += ll;
        }
      }
    }
    return cnt;
  }

  public void testCount() {
    assertEquals(120, lCount(AlleleSetProbabilityDiploid.LOOKUP));
  }
}
