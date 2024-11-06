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

import com.rtg.variant.bayes.Description;
import com.rtg.variant.util.arithmetic.SimplePossibility;

import junit.framework.TestCase;

/**
 */
public class HypothesesPriorTest extends TestCase {

  public void testHaploid() {
    final Description descr = new DescriptionCommon("X", "Y", "Z");
    final HypothesesPrior<Description> hyp = new HypothesesPrior<>(descr, SimplePossibility.SINGLETON, new double[] {0.1, 0.2, 0.7}, true, 2);
    //System.err.println(IntegralAbstract.toString(hyp));
    assertEquals(0.1, hyp.p(0), 0.001);
    assertEquals(0.2, hyp.p(1), 0.001);
    assertEquals(0.7, hyp.p(2), 0.001);
    assertEquals(2, hyp.reference());
  }

  public void testDiploid() {
    final double[] probs = {0.15, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.05};
    final HypothesesPrior<Description> hyp = new HypothesesPrior<>(DescriptionSnp.SINGLETON, SimplePossibility.SINGLETON, probs, true, 2);
    assertEquals(0.15, hyp.p(0), 0.001);
    assertEquals(0.10, hyp.p(1), 0.001);
    assertEquals(0.10, hyp.p(2), 0.001);
    assertEquals(0.10, hyp.p(3), 0.001);

    assertEquals(0.10,  hyp.p(4), 0.001); //A:C
    assertEquals(0.10,  hyp.p(5), 0.001); //C:G
    assertEquals(0.10, hyp.p(6), 0.001); //G:T
    assertEquals(0.10, hyp.p(7), 0.001); //A:G
    assertEquals(0.10,  hyp.p(8), 0.001); //C:T
    assertEquals(0.05, hyp.p(9), 0.001); //A:T
    assertEquals(2, hyp.reference());
  }

}
