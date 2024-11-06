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
