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

package com.rtg.variant.bayes;

import com.rtg.variant.bayes.snp.DescriptionCommon;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;
import com.rtg.variant.util.arithmetic.SimplePossibility;

import junit.framework.TestCase;

/**
 */
public class UnitFactorTest extends TestCase {

  public void test() {
    final Description desc = new DescriptionCommon("A", "B");
    final PossibilityArithmetic arith = SimplePossibility.SINGLETON;
    final Hypotheses<?> hyp = new MockHypotheses<>(desc, arith, false, new double[] {0.0, 0.0, 0.0}, 0);
    final Factor<?> factor = new UnitFactor<>(hyp, arith, hyp.size());
    assertTrue(hyp == factor.hypotheses());
    assertTrue(arith == factor.arithmetic());
    for (int i = 0; i < 3; ++i) {
      assertEquals(arith.one(), factor.p(i));
    }

  }
}
