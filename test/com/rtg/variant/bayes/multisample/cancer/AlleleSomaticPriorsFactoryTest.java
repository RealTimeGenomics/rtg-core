/*
 * Copyright (c) 2017. Real Time Genomics Limited.
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

import java.io.IOException;
import java.util.Arrays;

import com.rtg.util.InvalidParamsException;
import com.rtg.variant.bayes.Hypotheses;
import com.rtg.variant.bayes.HypothesesPowerSet;
import com.rtg.variant.bayes.MockHypotheses;
import com.rtg.variant.bayes.snp.DescriptionCommon;
import com.rtg.variant.util.arithmetic.SimplePossibility;

import junit.framework.TestCase;

/**
 */
public class AlleleSomaticPriorsFactoryTest extends TestCase {

  public void test() throws InvalidParamsException, IOException {
    final DescriptionCommon desc = new DescriptionCommon("", "A", "B");
    final Hypotheses<DescriptionCommon> hyp = new MockHypotheses<>(desc, SimplePossibility.SINGLETON, true, new double[] {0.0, 0.0, 0.0}, 0);
    final double[][] qf = new AlleleSomaticPriorsFactory<>(hyp).somaticQ(0.5);
    final double[][] q = SomaticPriorsAllele.makeQ(0.5, hyp, new HypothesesPowerSet<>(desc, SimplePossibility.SINGLETON, 0));
    assertTrue(Arrays.deepEquals(q, qf));
  }
}
