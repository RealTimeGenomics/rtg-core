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

package com.rtg.variant.bayes.multisample.cancer;

import com.rtg.variant.GenomePriorParams;
import com.rtg.variant.bayes.snp.HypothesesSnp;
import com.rtg.variant.util.arithmetic.SimplePossibility;

import junit.framework.TestCase;

/**
 */
public class HypothesesCancerTest extends TestCase {

  public void testHaploid() {
    final GenomePriorParams params = GenomePriorParams.builder().create();
    final HypothesesSnp snp = new HypothesesSnp(SimplePossibility.SINGLETON, params, true, 1);
    final HypothesesCancer hyp = new HypothesesCancer(snp, SimplePossibility.SINGLETON);
    assertFalse(hyp.haploid());
    assertEquals(1, hyp.reference());
    assertEquals(snp, hyp.subHypotheses());
  }

  public void testDiploid() {
    final GenomePriorParams params = GenomePriorParams.builder().create();
    final HypothesesSnp snp = new HypothesesSnp(SimplePossibility.SINGLETON, params, false, 2);
    final HypothesesCancer hyp = new HypothesesCancer(snp, SimplePossibility.SINGLETON);
    assertFalse(hyp.haploid());
    assertEquals(2, hyp.reference());
    assertEquals(snp, hyp.subHypotheses());
  }
}
