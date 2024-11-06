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

import com.rtg.util.integrity.Exam;
import com.rtg.variant.bayes.Code;
import com.rtg.variant.bayes.CodeDiploid;

import junit.framework.TestCase;

/**
 */
public class MendelianAlleleProbabilityHDDTest extends TestCase {

  private static void checkChildren(final double[][][][] lookup) {
    for (int j = 0; j < lookup.length; ++j) {
      for (int k = 0; k < lookup[j].length; ++k) {
        double sum = 0.0;
        for (int l = 0; l < lookup[j][k].length; ++l) {
          for (int m = 0; m < lookup[j][k][l].length; ++m) {
            if (l <= m) {
              sum += Math.exp(lookup[j][k][l][m]);
            }
            if (l != m) {
              Exam.assertEquals(lookup[j][k][l][m], lookup[j][k][m][l]);
            }
          }
        }
        Exam.assertEquals("j=" + j + " k=" + k, 1.0, sum, 0.0000001);
      }
    }
  }

  public void testLookup() {
    checkChildren(MendelianAlleleProbabilityHDD.LOOKUP);
  }

  public void testProbabilityLnLookup() {
    final Code code = new CodeDiploid(4);
    assertEquals(0.0, MendelianAlleleProbabilityHDD.SINGLETON_HD.probabilityLn(code, code.code(0), code.code(0, 0), code.code(0, 0)));

    assertEquals(Math.log(0.5), MendelianAlleleProbabilityHDD.SINGLETON_HD.probabilityLn(code, code.code(0), code.code(0, 1), code.code(0, 0)));
    assertEquals(Math.log(0.5), MendelianAlleleProbabilityHDD.SINGLETON_HD.probabilityLn(code, code.code(0), code.code(0, 1), code.code(0, 1)));

    assertEquals(Math.log(0.5), MendelianAlleleProbabilityHDD.SINGLETON_HD.probabilityLn(code, code.code(0), code.code(1, 2), code.code(0, 1)));
    assertEquals(Math.log(0.5), MendelianAlleleProbabilityHDD.SINGLETON_HD.probabilityLn(code, code.code(0), code.code(1, 2), code.code(0, 2)));

    assertEquals(Double.NEGATIVE_INFINITY, MendelianAlleleProbabilityHDD.SINGLETON_HD.probabilityLn(code, code.code(0), code.code(2, 3), code.code(0, 0)));

    assertEquals(Double.NEGATIVE_INFINITY, MendelianAlleleProbabilityHDD.SINGLETON_HD.probabilityLn(code, code.code(0), code.code(2, 3), code.code(1, 2)));
    assertEquals(Double.NEGATIVE_INFINITY, MendelianAlleleProbabilityHDD.SINGLETON_HD.probabilityLn(code, code.code(0), code.code(2, 3), code.code(3, 2)));

    assertEquals(Double.NEGATIVE_INFINITY, MendelianAlleleProbabilityHDD.SINGLETON_HD.probabilityLn(code, code.code(0), code.code(1, 2), code.code(3, 3)));
  }

  //test symmetries and the two ways of accessing the table
  public void testLookupCode() {
    final Code code = new CodeDiploid(4);
    for (int i = 0; i < code.rangeSize(); ++i)  {
      for (int j = 0; j < code.size(); ++j)  {
        //System.err.println("i=" + i + " j=" + j);
        for (int k = 0; k < code.size(); ++k)  {
          final double exp = MendelianAlleleProbabilityHDD.SINGLETON_HD.probabilityLn(code, i, j, k);
          assertEquals(exp, MendelianAlleleProbabilityHDD.SINGLETON_DH.probabilityLn(code, j, i, k));
          if (code.homozygous(j)) {
            assertEquals(exp, MendelianAlleleProbabilityHDD.SINGLETON_HD.probabilityLn(code, j, i, k));
            assertEquals(exp, MendelianAlleleProbabilityHDD.SINGLETON_DH.probabilityLn(code, i, j, k));
          }
        }
      }
    }
  }
}
