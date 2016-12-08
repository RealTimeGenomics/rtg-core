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
public class AlleleSetProbabilityHaploidTest extends TestCase {

  public void testLookup() {
    final Code code = new CodeDiploid(4);
    check(1.0 * 1.0 / 1.0, code, 0, code.code(0), code.code(0, 0));

    check(2.0 * 2.0 / 14.0, code, 0, code.code(0), code.code(0, 1));
    check(1.0 * 1.0 / 14.0, code, 1, code.code(0), code.code(0, 0));

    check(3.0 * 2.0 / 36.0, code, 0, code.code(0), code.code(1, 2));
    check(1.0 * 1.0 / 36.0, code, 2, code.code(0), code.code(1, 1));

    check(1.0 * 2.0 / 24.0, code, 3, code.code(0), code.code(1, 2));
  }

  private void check(final double exp, final Code code, final int ref, final int f, final int m) {
    assertEquals(Math.log(exp), AlleleSetProbabilityHaploid.SINGLETON_HD.probabilityLn(code, ref, f, m), 1e-7);
    assertEquals(Math.log(exp), AlleleSetProbabilityHaploid.SINGLETON_DH.probabilityLn(code, ref, m, f), 1e-7);
    if (code.homozygous(m)) {
      assertEquals(Math.log(exp), AlleleSetProbabilityHaploid.SINGLETON_DH.probabilityLn(code, ref, f, m), 1e-7);
      assertEquals(Math.log(exp), AlleleSetProbabilityHaploid.SINGLETON_HD.probabilityLn(code, ref, m, f), 1e-7);
    }
  }

  //test symmetries and the two ways of accessing the table
  public void testLookupCode() {
    final Code code = new CodeDiploid(4);
    for (int i = 0; i < code.rangeSize(); ++i)  {
      for (int j = 0; j < code.size(); ++j)  {
        //System.err.println("i=" + i + " j=" + j);
        final byte ref = 1;
        final double exp = AlleleSetProbabilityHaploid.SINGLETON_HD.probabilityLn(code, ref, i, j);
        assertEquals(exp, AlleleSetProbabilityHaploid.SINGLETON_DH.probabilityLn(code, ref, j, i));
        if (code.homozygous(j)) {
          assertEquals(exp, AlleleSetProbabilityHaploid.SINGLETON_HD.probabilityLn(code, ref, j, i));
          assertEquals(exp, AlleleSetProbabilityHaploid.SINGLETON_DH.probabilityLn(code, ref, i, j));
        }
      }
    }
  }

  public void testG() {
    assertEquals(Math.pow(4, 4), AlleleSetProbabilityHaploid.gSum());
    assertEquals(4, AlleleSetProbabilityHaploid.G.length);
  }


  static int lCount(final double[][][] lookup) {
    int cnt = 0;
    final int il0 = lookup[0].length;
    for (double[][] aLookup : lookup) {
      final int il = aLookup.length;
      assertEquals(il0, il);
      final int kl0 = aLookup[0].length;
      for (final double[] anALookup : aLookup) {
        final int kl = anALookup.length;
        assertEquals(kl0, kl);
        cnt += kl;
      }
    }
    return cnt;
  }

  public void testCount() {
    assertEquals(24, lCount(AlleleSetProbabilityHaploid.LOOKUP));
  }

}
