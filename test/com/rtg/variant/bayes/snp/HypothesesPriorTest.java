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

package com.rtg.variant.bayes.snp;

import com.rtg.variant.bayes.Description;
import com.rtg.variant.util.arithmetic.SimplePossibility;

import junit.framework.TestCase;

/**
 */
public class HypothesesPriorTest extends TestCase {

  public void testHaploid() {
    final Description descr = new DescriptionCommon("X", "Y", "Z");
    final HypothesesPrior<Description> hyp = new HypothesesPrior<>(descr, SimplePossibility.SINGLETON, new double[] {0.1, 0.2, 0.7}, true, 2);
    //System.err.println(IntegralAbstract.toString(hyp));
    assertEquals(0.1, hyp.p(0), 0.001);
    assertEquals(0.2, hyp.p(1), 0.001);
    assertEquals(0.7, hyp.p(2), 0.001);
    assertEquals(2, hyp.reference());
  }

  public void testDiploid() {
    final double[] probs = {0.15, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.05};
    final HypothesesPrior<Description> hyp = new HypothesesPrior<>(DescriptionSnp.SINGLETON, SimplePossibility.SINGLETON, probs, true, 2);
    assertEquals(0.15, hyp.p(0), 0.001);
    assertEquals(0.10, hyp.p(1), 0.001);
    assertEquals(0.10, hyp.p(2), 0.001);
    assertEquals(0.10, hyp.p(3), 0.001);

    assertEquals(0.10,  hyp.p(4), 0.001); //A:C
    assertEquals(0.10,  hyp.p(5), 0.001); //C:G
    assertEquals(0.10, hyp.p(6), 0.001); //G:T
    assertEquals(0.10, hyp.p(7), 0.001); //A:G
    assertEquals(0.10,  hyp.p(8), 0.001); //C:T
    assertEquals(0.05, hyp.p(9), 0.001); //A:T
    assertEquals(2, hyp.reference());
  }

}
