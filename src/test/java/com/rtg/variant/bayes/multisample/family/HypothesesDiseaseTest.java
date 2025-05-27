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

import com.rtg.variant.bayes.snp.DescriptionSnp;
import com.rtg.variant.util.arithmetic.LogApproximatePossibility;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;

import junit.framework.TestCase;

/**
 */
public class HypothesesDiseaseTest extends TestCase {

  public void test() {
    final PossibilityArithmetic arith = LogApproximatePossibility.SINGLETON;
    final HypothesesDisease hyp = new HypothesesDisease(DescriptionSnp.SINGLETON, 0.1, 3);

    assertEquals(0.1, arith.poss2Prob(hyp.prior(0)), 0.000001);
    final double dis = (1.0 - 0.1) / 3;
    assertEquals(dis, arith.poss2Prob(hyp.prior(1)), 0.00001);
    assertEquals(dis, arith.poss2Prob(hyp.prior(2)), 0.00001);
    assertEquals(dis, arith.poss2Prob(hyp.prior(3)), 0.00001);
    assertEquals(0.0, arith.poss2Prob(hyp.prior(4)), 0.00001);

    assertEquals(4, hyp.reference());
    assertEquals(5, hyp.size());
  }
}
