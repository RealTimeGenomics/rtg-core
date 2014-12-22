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

import com.rtg.util.integrity.Exam;
import com.rtg.variant.bayes.Code;
import com.rtg.variant.bayes.CodeDiploid;

import junit.framework.TestCase;

/**
 */
public class MendelianAlleleProbabilityHDHTest extends TestCase {

  private static void checkChildren(final double[][] lookup) {
    for (int i = 0; i < lookup.length; i++) {
      double sum = 0.0;
      for (int j = 0; j < lookup[i].length; j++) {
        sum += Math.exp(lookup[i][j]);
      }
      Exam.assertEquals("j=" + i, 1.0, sum, 0.0000001);
    }
  }

  public void testLookup() {
    checkChildren(MendelianAlleleProbabilityHDH.LOOKUP);
  }

  public void testProbabilityLnLookup() {
    final Code code = new CodeDiploid(4);
    check(1.0, code, code.code(0), code.code(0, 0), code.code(0));

    check(0.5, code, code.code(0), code.code(0, 1), code.code(0));
    check(0.5, code, code.code(0), code.code(0, 1), code.code(1));

    check(0.5, code, code.code(0), code.code(1, 2), code.code(1));
    check(0.5, code, code.code(0), code.code(1, 2), code.code(2));

    check(1.0, code, code.code(0), code.code(1, 1), code.code(1));

    check(0.0, code, code.code(0), code.code(2, 3), code.code(0, 0));

    check(0.0, code, code.code(0), code.code(2, 3), code.code(1));
    check(0.0, code, code.code(0), code.code(2, 3), code.code(0));

    check(0.0, code, code.code(0), code.code(1, 2), code.code(3));
  }

  private void check(final double exp, final Code code, final int f, final int m, final int c) {
    assertEquals(Math.log(exp), MendelianAlleleProbabilityHDH.SINGLETON_HD.probabilityLn(code, f, m, c));
  }

  //test symmetries and the two ways of accessing the table
  public void testLookupCode() {
    final Code code = new CodeDiploid(4);
    for (int i = 0; i < code.rangeSize(); i++)  {
      for (int j = 0; j < code.size(); j++)  {
        //System.err.println("i=" + i + " j=" + j);
        for (int k = 0; k < code.rangeSize(); k++)  {
          final double exp = MendelianAlleleProbabilityHDH.SINGLETON_HD.probabilityLn(code, i, j, k);
          assertEquals(exp, MendelianAlleleProbabilityHDH.SINGLETON_DH.probabilityLn(code, j, i, k));
        }
      }
    }
  }
}
