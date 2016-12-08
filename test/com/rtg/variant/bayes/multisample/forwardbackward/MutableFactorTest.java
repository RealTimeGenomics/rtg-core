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
