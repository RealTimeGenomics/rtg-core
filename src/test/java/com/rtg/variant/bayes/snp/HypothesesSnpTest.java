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

package com.rtg.variant.bayes.snp;

import com.rtg.variant.GenomePriorParams;
import com.rtg.variant.util.arithmetic.LogPossibility;

import junit.framework.TestCase;

/**
 */
public class HypothesesSnpTest extends TestCase {

  //Taken from HaploidSnpBayesianTest.testNewPriors()
  public void testHaploid() {
    final GenomePriorParams params = GenomePriorParams.builder().create();
    final HypothesesSnp hyp = new HypothesesSnp(LogPossibility.SINGLETON, params, true, 1);
    //System.err.println(IntegralAbstract.toString(hyp));
    assertEquals(-9.257, hyp.p(0), 0.001);
    assertEquals(-0.000, hyp.p(1), 0.001);
    assertEquals(-9.172, hyp.p(2), 0.001);
    assertEquals(-7.884, hyp.p(3), 0.001);
    assertEquals(1, hyp.reference());
  }

  //Taken from HaploidSnpBayesianTest.testNewPriors()
  public void testDiploid() {
    final GenomePriorParams params = GenomePriorParams.builder().create();
    final HypothesesSnp hyp = new HypothesesSnp(LogPossibility.SINGLETON, params, false, 1);
    //    System.err.println(IntegralAbstract.toString(hyp));
    //    for (int i = 0; i < hyp.size(); ++i) {
    //      System.err.println("i=" + i + " name=" + hyp.name(i));
    //    }
    assertEquals(-9.257, hyp.p(0), 0.001);
    assertEquals(-0.002, hyp.p(1), 0.001);
    assertEquals(-9.172, hyp.p(2), 0.001);
    assertEquals(-7.885, hyp.p(3), 0.001);

    assertEquals(-8.724,  hyp.p(4), 0.001); //A:C
    assertEquals(-8.724,  hyp.p(5), 0.001); //C:G
    assertEquals(-14.894, hyp.p(6), 0.001); //G:T
    assertEquals(-14.894, hyp.p(7), 0.001); //A:G
    assertEquals(-7.288,  hyp.p(8), 0.001); //C:T
    assertEquals(-14.894, hyp.p(9), 0.001); //A:T
    assertEquals(1, hyp.reference());
  }
}
