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

import com.rtg.variant.GenomePriorParams;
import com.rtg.variant.util.arithmetic.LogPossibility;

import junit.framework.TestCase;

/**
 */
public class HypothesesSnpTest extends TestCase {

  //Taken from HaploidSnpBayesianTest.testNewPriors()
  public void testHaploid() {
    final GenomePriorParams params = GenomePriorParams.builder().create();
    final HypothesesSnp hyp = new HypothesesSnp(LogPossibility.SINGLETON, params, true, 1);
    //System.err.println(IntegralAbstract.toString(hyp));
    assertEquals(-9.257, hyp.p(0), 0.001);
    assertEquals(-0.000, hyp.p(1), 0.001);
    assertEquals(-9.172, hyp.p(2), 0.001);
    assertEquals(-7.884, hyp.p(3), 0.001);
    assertEquals(1, hyp.reference());
  }

  //Taken from HaploidSnpBayesianTest.testNewPriors()
  public void testDiploid() {
    final GenomePriorParams params = GenomePriorParams.builder().create();
    final HypothesesSnp hyp = new HypothesesSnp(LogPossibility.SINGLETON, params, false, 1);
    //    System.err.println(IntegralAbstract.toString(hyp));
    //    for (int i = 0; i < hyp.size(); ++i) {
    //      System.err.println("i=" + i + " name=" + hyp.name(i));
    //    }
    assertEquals(-9.257, hyp.p(0), 0.001);
    assertEquals(-0.002, hyp.p(1), 0.001);
    assertEquals(-9.172, hyp.p(2), 0.001);
    assertEquals(-7.885, hyp.p(3), 0.001);

    assertEquals(-8.724,  hyp.p(4), 0.001); //A:C
    assertEquals(-8.724,  hyp.p(5), 0.001); //C:G
    assertEquals(-14.894, hyp.p(6), 0.001); //G:T
    assertEquals(-14.894, hyp.p(7), 0.001); //A:G
    assertEquals(-7.288,  hyp.p(8), 0.001); //C:T
    assertEquals(-14.894, hyp.p(9), 0.001); //A:T
    assertEquals(1, hyp.reference());
  }
}
