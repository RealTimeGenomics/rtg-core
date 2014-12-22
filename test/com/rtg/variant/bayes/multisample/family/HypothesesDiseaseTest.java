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
