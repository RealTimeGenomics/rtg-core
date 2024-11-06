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
public class MendelianAlleleProbabilityHDHTest extends TestCase {

  private static void checkChildren(final double[][] lookup) {
    for (int i = 0; i < lookup.length; ++i) {
      double sum = 0.0;
      for (int j = 0; j < lookup[i].length; ++j) {
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
    for (int i = 0; i < code.rangeSize(); ++i)  {
      for (int j = 0; j < code.size(); ++j)  {
        //System.err.println("i=" + i + " j=" + j);
        for (int k = 0; k < code.rangeSize(); ++k)  {
          final double exp = MendelianAlleleProbabilityHDH.SINGLETON_HD.probabilityLn(code, i, j, k);
          assertEquals(exp, MendelianAlleleProbabilityHDH.SINGLETON_DH.probabilityLn(code, j, i, k));
        }
      }
    }
  }
}
