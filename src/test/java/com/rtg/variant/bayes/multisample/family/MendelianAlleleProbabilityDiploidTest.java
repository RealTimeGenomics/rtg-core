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
public class MendelianAlleleProbabilityDiploidTest extends TestCase {

  private static void checkChildren(final double[][][][][] lookup) {
    for (int i = 0; i < lookup.length; ++i) {
      for (int j = 0; j < lookup[i].length; ++j) {
        for (int k = 0; k < lookup[i][j].length; ++k) {
          double sum = 0.0;
          for (int l = 0; l < lookup[i][j][k].length; ++l) {
            for (int m = 0; m < lookup[i][j][k][l].length; ++m) {
              if (l <= m) {
                sum += Math.exp(lookup[i][j][k][l][m]);
              }
              if (l != m) {
                Exam.assertEquals(lookup[i][j][k][l][m], lookup[i][j][k][m][l]);
              }
            }
          }
          Exam.assertEquals("i=" + i + " j=" + j + " k=" + k, 1.0, sum, 0.0000001);
        }
      }
    }
  }

  public void testLookup() {
    checkChildren(MendelianAlleleProbabilityDiploid.LOOKUP);
  }

  public void testProbabilityLnLookup() {
    final Code code = new CodeDiploid(4);
    assertEquals(0.0, MendelianAlleleProbabilityDiploid.SINGLETON.probabilityLn(code, code.code(0), code.code(1), code.code(0, 1)));
    assertEquals(Math.log(0.5), MendelianAlleleProbabilityDiploid.SINGLETON.probabilityLn(code, code.code(0, 1), code.code(1), code.code(0, 1)));
    assertEquals(Math.log(0.25), MendelianAlleleProbabilityDiploid.SINGLETON.probabilityLn(code, code.code(0, 1), code.code(0, 1), code.code(0)));
    assertEquals(Math.log(0.25), MendelianAlleleProbabilityDiploid.SINGLETON.probabilityLn(code, code.code(0, 1), code.code(2, 3), code.code(0, 3)));
    assertEquals(Double.NEGATIVE_INFINITY, MendelianAlleleProbabilityDiploid.SINGLETON.probabilityLn(code, code.code(0, 1), code.code(2, 3), code.code(0)));

    assertEquals(Double.NEGATIVE_INFINITY, MendelianAlleleProbabilityDiploid.SINGLETON.probabilityLn(code, code.code(0), code.code(1), code.code(2, 0)));
    assertEquals(Double.NEGATIVE_INFINITY, MendelianAlleleProbabilityDiploid.SINGLETON.probabilityLn(code, code.code(0), code.code(2), code.code(2, 1)));
  }

  //test symmetries and the two ways of accessing the table
  public void testLookupCode() {
    final Code code = new CodeDiploid(4);
    for (int i = 0; i < code.size(); ++i)  {
      for (int j = 0; j < code.size(); ++j)  {
        //System.err.println("i=" + i + " j=" + j);
        for (int k = 0; k < code.size(); ++k)  {
          final double exp = MendelianAlleleProbabilityDiploid.SINGLETON.probabilityLn(code, i, j, k);
          assertEquals(exp, MendelianAlleleProbabilityDiploid.SINGLETON.probabilityLn(code, j, i, k));
        }
      }
    }
  }
}
