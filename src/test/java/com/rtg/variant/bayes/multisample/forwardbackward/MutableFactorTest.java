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
package com.rtg.variant.bayes.multisample.forwardbackward;

import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.Hypotheses;
import com.rtg.variant.bayes.MockHypotheses;
import com.rtg.variant.bayes.snp.DescriptionCommon;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;
import com.rtg.variant.util.arithmetic.SimplePossibility;

import junit.framework.TestCase;

/**
 */
public class MutableFactorTest extends TestCase {

  public void test() {
    final Description desc = new DescriptionCommon("A", "B");
    final PossibilityArithmetic arith = SimplePossibility.SINGLETON;
    final Hypotheses<?> hyp = new MockHypotheses<>(desc, arith, false, new double[] {0.0, 0.3, 0.5}, 0);
    final MutableFactor<?> mhv = new MutableFactor<>(hyp, arith, hyp.size());
    assertTrue(mhv.hypotheses() == hyp);
    assertTrue(mhv.arithmetic() == arith);
    for (int i = 0; i < 3; ++i) {
      assertEquals(arith.zero(), mhv.p(i));
    }

    mhv.set(0, arith.prob2Poss(0.1));
    mhv.set(2, arith.prob2Poss(0.3));
    final double[] values = {0.1, 0.0, 0.3};
    for (int i = 0; i < 3; ++i) {
      assertEquals(values[i], arith.poss2Prob(mhv.p(i)));
    }
    try {
      mhv.p(3);
      fail();
    } catch (final RuntimeException e) {
      //expected
    }
    try {
      mhv.p(-1);
      fail();
    } catch (final RuntimeException e) {
      //expected
    }

    mhv.set(1, arith.prob2Poss(0.2));
    assertEquals(0.2, arith.poss2Prob(mhv.p(1)));
  }

}
