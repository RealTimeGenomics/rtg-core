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
import com.rtg.variant.bayes.snp.StatisticsSnp;
import com.rtg.variant.util.arithmetic.SimplePossibility;

import junit.framework.TestCase;

/**
 * @author kurt
 */
public class NoAlleleBalanceTest extends TestCase {
  public void test() {
    final DescriptionCommon descriptionCommon = new DescriptionCommon("A", "C");
    assertEquals(0.0, new NoAlleleBalance().alleleBalanceLn(0, new MockHypotheses<Description>(descriptionCommon, SimplePossibility.SINGLETON, true, new double[] {0.1, 0.1}, 0), new StatisticsSnp(descriptionCommon)));
  }

}