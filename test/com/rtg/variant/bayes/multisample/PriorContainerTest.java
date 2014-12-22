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
package com.rtg.variant.bayes.multisample;


import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.MockHypotheses;
import com.rtg.variant.bayes.UnitFactor;
import com.rtg.variant.bayes.multisample.forwardbackward.BContainer;
import com.rtg.variant.bayes.snp.DescriptionCommon;
import com.rtg.variant.bayes.snp.HypothesesNone;
import com.rtg.variant.bayes.snp.HypothesesPrior;
import com.rtg.variant.util.arithmetic.SimplePossibility;

import junit.framework.TestCase;

/**
 */
public class PriorContainerTest extends TestCase {

  public void test() {
    final HypothesesPrior<Description> hap = new MockHypotheses<Description>(new DescriptionCommon("a"), SimplePossibility.SINGLETON, true, new double[] {1.0}, 0);
    final HypothesesPrior<Description> dip = new MockHypotheses<Description>(new DescriptionCommon("a"), SimplePossibility.SINGLETON, false, new double[] {1.0}, 0);
    final BContainer[] bs = {new BContainer(new UnitFactor<>(hap, hap.arithmetic(), 1))};
    final HaploidDiploidHypotheses<HypothesesPrior<Description>> hdp = new HaploidDiploidHypotheses<>(HypothesesNone.SINGLETON, hap, dip, false, null);
    final PriorContainer<?> pc = new PriorContainer<>(hdp, bs);
    assertEquals(hdp, pc.getHypotheses());
    assertTrue(bs == pc.getBs());
  }
}
