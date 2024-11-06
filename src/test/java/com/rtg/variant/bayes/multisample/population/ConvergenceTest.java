/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
    for (int i = 0; i < 100; ++i) {
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
    for (int i = 0; i < 100; ++i) {
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
    for (int i = 0; i < 100; ++i) {
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
    for (int i = 0; i < 100; ++i) {
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
