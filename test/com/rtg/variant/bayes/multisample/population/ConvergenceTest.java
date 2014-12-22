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

package com.rtg.variant.bayes.multisample.population;


import java.util.Arrays;
import java.util.Random;

import com.rtg.variant.bayes.CodeDiploid;
import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.Hypotheses;
import com.rtg.variant.bayes.MockHypotheses;
import com.rtg.variant.bayes.snp.DescriptionSnp;
import com.rtg.variant.bayes.snp.HypothesesPrior;
import com.rtg.variant.util.arithmetic.SimplePossibility;

import junit.framework.TestCase;

/**
 */
public class ConvergenceTest extends TestCase {

  public void testError() {
    final Description descr = DescriptionSnp.SINGLETON;
    final Hypotheses<?> hyp = new MockHypotheses<>(descr, SimplePossibility.SINGLETON, true, new double[] {0.25, 0.25, 0.25, 0.25}, 0);
    final Random random = new Random(43);
    final int[] counts = new int[4];
    for (int i = 0; i < 100; i++) {
      final int r = Convergence.error(0.4, 1, hyp, random);
      //System.err.println("r=" + r);
      counts[r]++;
    }
    final String str = Arrays.toString(counts);
    //System.err.println(str);
    assertEquals("[13, 59, 13, 15]", str);
  }

  public void testChoose() {
    final Description descr = DescriptionSnp.SINGLETON;
    final Hypotheses<?> hyp = new MockHypotheses<>(descr, SimplePossibility.SINGLETON, false, AbstractEstimatorTest.uniform(10), 0);
    final Random random = new Random(43);
    final int[] counts = new int[10];
    for (int i = 0; i < 100; i++) {
      final int r = Convergence.choose(5, hyp, random);
      //System.err.println("r=" + r);
      counts[r]++;
    }
    final String str = Arrays.toString(counts);
    //System.err.println(str);
    assertEquals("[0, 50, 50, 0, 0, 0, 0, 0, 0, 0]", str);
  }

  public void testSampleHaploid() {
    final double[] priors = {0.1, 0.2, 0.3, 0.4};
    final Random random = new Random(43);
    final int[] counts = new int[4];
    for (int i = 0; i < 100; i++) {
      final int r = Convergence.sampleHaploid(priors, random);
      //System.err.println("r=" + r);
      counts[r]++;
    }
    final String str = Arrays.toString(counts);
    //System.err.println(str);
    assertEquals("[8, 23, 31, 38]", str);
  }

  public void testSample() {
    final double[] priors = {0.1, 0.2, 0.3, 0.4};
    final Random random = new Random(43);
    final CodeDiploid code = new CodeDiploid(4);
    final int[] counts = new int[10];
    for (int i = 0; i < 100; i++) {
      final int r = Convergence.sample(priors, code, random);
      //System.err.println("r=" + r);
      counts[r]++;
    }
    final String str = Arrays.toString(counts);
    //System.err.println(str);
    assertEquals("[1, 3, 10, 16, 4, 11, 28, 4, 16, 7]", str);
  }

  public void testDistr() {
    final Description descr = DescriptionSnp.SINGLETON;
    final HypothesesPrior<?> hahyp = new MockHypotheses<>(descr, SimplePossibility.SINGLETON, true, new double[] {0.1, 0.2, 0.3, 0.4}, 0);
    final double[] distr = Convergence.distr(hahyp);
    assertEquals(4, distr.length);
    assertEquals(0.1, distr[0]);
    assertEquals(0.2, distr[1]);
    assertEquals(0.3, distr[2]);
    assertEquals(0.4, distr[3]);
  }
}
