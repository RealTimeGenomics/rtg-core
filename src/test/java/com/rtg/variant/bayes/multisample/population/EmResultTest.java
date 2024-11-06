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

package com.rtg.variant.bayes.multisample.population;


import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.MockGenotypeMeasure;
import com.rtg.variant.bayes.MockHypotheses;
import com.rtg.variant.bayes.multisample.HaploidDiploidHypotheses;
import com.rtg.variant.bayes.multisample.HypothesisScore;
import com.rtg.variant.bayes.multisample.HypothesisScores;
import com.rtg.variant.bayes.snp.DescriptionCommon;
import com.rtg.variant.bayes.snp.HypothesesNone;
import com.rtg.variant.bayes.snp.HypothesesPrior;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;
import com.rtg.variant.util.arithmetic.SimplePossibility;

import junit.framework.TestCase;

/**
 */
public class EmResultTest extends TestCase {

  public void test() {
    final PossibilityArithmetic arith = SimplePossibility.SINGLETON;
    final DescriptionCommon descr = new DescriptionCommon("X", "Y", "Z");
    final MockHypotheses<Description> haploid = new MockHypotheses<>(descr, arith, true, AbstractEstimatorTest.uniform(3), 0);
    final MockHypotheses<Description> diploid = new MockHypotheses<>(descr, arith, false, AbstractEstimatorTest.uniform(6), 0);
    final HypothesisScore[] calls = new HypothesisScore[42];
    final HaploidDiploidHypotheses<HypothesesPrior<Description>> hdh = new HaploidDiploidHypotheses<>(HypothesesNone.SINGLETON, haploid, diploid);
    final EmResult<HypothesesPrior<Description>> emr = new EmResult<>(hdh, new HypothesisScores(calls, false, 0, null), null);
    assertTrue(emr.calls().getScores() == calls);
    assertTrue(emr.getPriorContainer().getHypotheses().haploid() == haploid);
    assertTrue(emr.getPriorContainer().getHypotheses().diploid() == diploid);
  }

  HypothesisScore[] scores(final int...hyp) {
    final HypothesisScore[] calls = new HypothesisScore[hyp.length];
    for (int i = 0; i < hyp.length; ++i) {
      calls[i] = new HypothesisScore(new MockGenotypeMeasure(0, hyp[i], 0.0, 0));
    }
    return calls;
  }

  public void testDifference() {
    final HypothesisScore[] calls0 = scores(1, 2, 3);
    final HypothesisScore[] calls1 = scores(4, 2, 3);
    final HypothesisScore[] calls2 = scores(1, 2, 4);
    final HypothesisScore[] calls3 = scores(4, 5, 6);
    assertEquals(0, EmResult.difference(calls0, calls0));
    assertEquals(1, EmResult.difference(calls0, calls1));
    assertEquals(1, EmResult.difference(calls0, calls2));
    assertEquals(3, EmResult.difference(calls0, calls3));
  }

  public void testDifference0() {
    final HypothesisScore[] calls0 = scores();
    assertEquals(0, EmResult.difference(calls0, calls0));
  }

  public void testDifference1() {
    final PossibilityArithmetic arith = SimplePossibility.SINGLETON;
    final DescriptionCommon descr = new DescriptionCommon("X", "Y", "Z");
    final MockHypotheses<Description> haploid = new MockHypotheses<>(descr, arith, true, AbstractEstimatorTest.uniform(3), 0);
    final MockHypotheses<Description> diploid = new MockHypotheses<>(descr, arith, false, AbstractEstimatorTest.uniform(6), 0);
    final HypothesisScores calls = new HypothesisScores(scores(1, 2, 3), false, 0, null);
    final HaploidDiploidHypotheses<HypothesesPrior<Description>> hdh = new HaploidDiploidHypotheses<>(HypothesesNone.SINGLETON, haploid, diploid);
    final EmResult<HypothesesPrior<Description>> emr = new EmResult<>(hdh, calls, null);
    final HypothesisScores calls1 = new HypothesisScores(scores(4, 2, 3), false, 0, null);
    assertEquals(1, emr.difference(calls1));
  }
}
