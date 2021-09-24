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
import com.rtg.variant.bayes.multisample.population.DescriptionCounts;
import com.rtg.variant.bayes.snp.DescriptionNone;
import com.rtg.variant.bayes.snp.HypothesesNone;
import com.rtg.variant.bayes.snp.HypothesesPrior;
import com.rtg.variant.util.arithmetic.SimplePossibility;

import junit.framework.TestCase;

/**
 */
public class HaploidDiploidHypothesesTest extends TestCase {

  public void testInitialization() {
    final HypothesesPrior<Description> hap = new HypothesesPrior<>(DescriptionNone.SINGLETON, SimplePossibility.SINGLETON, new double[]{}, true, 0);
    final HypothesesPrior<Description> dip = new HypothesesPrior<>(DescriptionNone.SINGLETON, SimplePossibility.SINGLETON, new double[]{}, false, 1);
    final DescriptionCounts dc = new DescriptionCounts(2, 0);
    final HaploidDiploidHypotheses<HypothesesPrior<Description>> hdh = new HaploidDiploidHypotheses<>(HypothesesNone.SINGLETON, hap, dip, true, dc);
    assertEquals(hap, hdh.haploid());
    assertEquals(dip, hdh.diploid());
    assertTrue(hdh.isDefault());
    assertNotNull(hdh.getDescriptionCounts());
  }
}
