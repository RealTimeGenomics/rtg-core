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
public class MendelianAlleleProbabilityHHDDeNovoTest extends TestCase {

  public void test() {
    final Code code = new CodeDiploid(4);
    assertEquals(Double.NEGATIVE_INFINITY, MendelianAlleleProbabilityHHDDeNovo.SINGLETON_HHD.probabilityLn(code, 0, 0, 0));
    assertEquals(Math.log(1.0 / 3), MendelianAlleleProbabilityHHDDeNovo.SINGLETON_HHD.probabilityLn(code, 0, 0, 1));
    assertEquals(Double.NEGATIVE_INFINITY, MendelianAlleleProbabilityHHDDeNovo.SINGLETON_HHD.probabilityLn(code, 1, 0, code.code(1, 0)));
    assertEquals(Math.log(1.0 / 3), MendelianAlleleProbabilityHHDDeNovo.SINGLETON_HHD.probabilityLn(code, 1, 1, code.code(1, 0)));
  }

  //test there is no overlap between non de novo and de novo tables
  public void testOverlap() {
    final Code code = new CodeDiploid(4);
    for (int j = 0; j < code.rangeSize(); ++j) {
      for (int k = 0; k < code.rangeSize(); ++k) {
        if (MendelianAlleleProbabilityHHD.SINGLETON_HHD.probabilityLn(code, j, -1, k) > Double.NEGATIVE_INFINITY) {
          assertEquals(Double.NEGATIVE_INFINITY, MendelianAlleleProbabilityHHDDeNovo.SINGLETON_HHD.probabilityLn(code, j, -1, k));
        }
        if (MendelianAlleleProbabilityHHD.SINGLETON_HHD.probabilityLn(code, -1, j, k) > Double.NEGATIVE_INFINITY) {
          assertEquals(Double.NEGATIVE_INFINITY, MendelianAlleleProbabilityHHDDeNovo.SINGLETON_HHD.probabilityLn(code, -1, j, k));
        }
      }
    }
  }
}
