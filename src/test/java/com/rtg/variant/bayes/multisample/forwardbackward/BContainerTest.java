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
import com.rtg.variant.bayes.MockHypotheses;
import com.rtg.variant.bayes.UnitFactor;
import com.rtg.variant.bayes.snp.DescriptionCommon;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;
import com.rtg.variant.util.arithmetic.SimplePossibility;

import junit.framework.TestCase;

/**
 */
public class BContainerTest extends TestCase {
  public void test() {
    final PossibilityArithmetic arith = SimplePossibility.SINGLETON;
    final BContainer b = new BContainer(new UnitFactor<>(new MockHypotheses<Description>(new DescriptionCommon("."), arith, true, new double[] {1.0}, 0), arith, 1));
    assertEquals(arith.one(), b.getB(0).p(0));
    assertEquals(1, b.size());
  }
}
