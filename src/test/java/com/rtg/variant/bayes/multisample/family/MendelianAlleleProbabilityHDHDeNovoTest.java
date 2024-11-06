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

import com.rtg.variant.bayes.Code;
import com.rtg.variant.bayes.CodeDiploid;

import junit.framework.TestCase;

/**
 */
public class MendelianAlleleProbabilityHDHDeNovoTest extends TestCase {

  public void testUniHypothesisLookup() {
    // Non de novo should give no likelihood
    final Code code = new CodeDiploid(1);
    assertEquals(Double.NEGATIVE_INFINITY, MendelianAlleleProbabilityHDHDeNovo.SINGLETON_HD.probabilityLn(code, code.code(0), code.code(0), code.code(0)));
  }

  public void testProbabilityLnLookup() {
    final Code code = new CodeDiploid(4);

    // Non de novo should give no likelihood
    assertEquals(Double.NEGATIVE_INFINITY, MendelianAlleleProbabilityHDHDeNovo.SINGLETON_HD.probabilityLn(code, code.code(0), code.code(1), code.code(1)));
    assertEquals(Double.NEGATIVE_INFINITY, MendelianAlleleProbabilityHDHDeNovo.SINGLETON_HD.probabilityLn(code, code.code(0), code.code(0, 1), code.code(1)));
    assertEquals(Double.NEGATIVE_INFINITY, MendelianAlleleProbabilityHDHDeNovo.SINGLETON_HD.probabilityLn(code, code.code(0), code.code(1, 2), code.code(2)));

    // De novo table entries
    assertEquals(Math.log(1.0 / 3), MendelianAlleleProbabilityHDHDeNovo.SINGLETON_HD.probabilityLn(code, code.code(0), code.code(0, 1), code.code(2)), 0.0001);
    assertEquals(Math.log(1.0 / 3), MendelianAlleleProbabilityHDHDeNovo.SINGLETON_HD.probabilityLn(code, code.code(0), code.code(1), code.code(2)), 0.0001);
  }

  //test symmetries and the two ways of accessing the table
  public void testLookupCode() {
    final Code code = new CodeDiploid(4);
    for (int i = 0; i < code.rangeSize(); ++i)  {
      for (int j = 0; j < code.size(); ++j)  {
        //System.err.println("i=" + i + " j=" + j);
        for (int k = 0; k < code.rangeSize(); ++k)  {
          final double exp = MendelianAlleleProbabilityHDHDeNovo.SINGLETON_HD.probabilityLn(code, i, j, k);
          assertEquals("f=" + code.a(i) + " m=" + code.a(j) + "/" + code.bc(j) + " c=" + code.a(k), exp, MendelianAlleleProbabilityHDHDeNovo.SINGLETON_DH.probabilityLn(code, j, i, k));
        }
      }
    }
  }

//test there is no overlap between non de novo and de novo tables
  public void testOverlap() {
    final Code code = new CodeDiploid(4);
    for (int i = 0; i < code.rangeSize(); ++i) {
      for (int j = 0; j < code.size(); ++j) {
        for (int k = 0; k < code.rangeSize(); ++k) {
          if (MendelianAlleleProbabilityHDH.SINGLETON_HD.probabilityLn(code, i, j, k) > Double.NEGATIVE_INFINITY) {
            assertEquals(Double.NEGATIVE_INFINITY, MendelianAlleleProbabilityHDHDeNovo.SINGLETON_HD.probabilityLn(code, i, j, k));
          }
          if (MendelianAlleleProbabilityHDH.SINGLETON_DH.probabilityLn(code, j, i, k) > Double.NEGATIVE_INFINITY) {
            assertEquals("father=" + code.a(i) + " mother=" + code.a(j) + "|" + code.bc(j) + " child=" + code.a(k) + "|" + code.bc(k), Double.NEGATIVE_INFINITY, MendelianAlleleProbabilityHDHDeNovo.SINGLETON_DH.probabilityLn(code, j, i, k));
          }
        }
      }
    }
  }
}
