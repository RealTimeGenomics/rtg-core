/*
 * Copyright (c) 2014. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
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
    final MockHypotheses<Description> haploid = new MockHypotheses<Description>(descr, arith, true, AbstractEstimatorTest.uniform(3), 0);
    final MockHypotheses<Description> diploid = new MockHypotheses<Description>(descr, arith, false, AbstractEstimatorTest.uniform(6), 0);
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
    final MockHypotheses<Description> haploid = new MockHypotheses<Description>(descr, arith, true, AbstractEstimatorTest.uniform(3), 0);
    final MockHypotheses<Description> diploid = new MockHypotheses<Description>(descr, arith, false, AbstractEstimatorTest.uniform(6), 0);
    final HypothesisScores calls = new HypothesisScores(scores(1, 2, 3), false, 0, null);
    final HaploidDiploidHypotheses<HypothesesPrior<Description>> hdh = new HaploidDiploidHypotheses<>(HypothesesNone.SINGLETON, haploid, diploid);
    final EmResult<HypothesesPrior<Description>> emr = new EmResult<>(hdh, calls, null);
    final HypothesisScores calls1 = new HypothesisScores(scores(4, 2, 3), false, 0, null);
    assertEquals(1, emr.difference(calls1));
  }
}
