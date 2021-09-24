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

import com.rtg.variant.GenomePriorParams;
import com.rtg.variant.bayes.snp.HypothesesSnp;
import com.rtg.variant.util.arithmetic.SimplePossibility;

import junit.framework.TestCase;

/**
 */
public class NoNonIdentityMeasureTest extends TestCase {
  public void test() {
    final NoNonIdentityMeasure measure = new NoNonIdentityMeasure(new ArrayGenotypeMeasure(SimplePossibility.SINGLETON, new double[] {0.4, 0.2, 0.1, 0.1, 0.2}, new HypothesesSnp(SimplePossibility.SINGLETON, GenomePriorParams.builder().create(), false, 2)));
    assertEquals(SimplePossibility.SINGLETON, measure.arithmetic());
    assertEquals(0, measure.best());
    assertEquals(2, measure.reference());
    assertEquals(5, measure.size());
    assertEquals(0.4 / 0.6, measure.bestPosterior());
    assertTrue(Double.isNaN(measure.nonIdentityPosterior()));
  }
}
