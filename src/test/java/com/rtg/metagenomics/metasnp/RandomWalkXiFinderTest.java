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
package com.rtg.metagenomics.metasnp;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.variant.util.arithmetic.LogPossibility;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;
import com.rtg.variant.util.arithmetic.SimplePossibility;

import junit.framework.TestCase;

/**
 */
public class RandomWalkXiFinderTest extends TestCase {
  
  @Override
  public void setUp() {
    Diagnostic.setLogStream();
  }

  private static final List<int[]> AT = Collections.singletonList(new int[] {0, 1});

  public void testEasyBalanced() {
    final RandomWalkXiFinder xiFinder = new RandomWalkXiFinder();
    final XiScore xis = xiFinder.maximizeSingleSample(AT, Collections.singletonList(new double[][] {{10, 10, 0, 0}}), 0);
    final double[] xi = xis.mXi;
    assertEquals(1, xi[0] + xi[1], 1e-10);
    assertEquals(0.5, xi[0], 0.1);
    assertEquals(0.5, xi[1], 0.1);
    assertEquals(9.536742780032061e-7, xis.mScore, 1e-17);
  }
  public void testEasyBalancedCG() {
    final RandomWalkXiFinder xiFinder = new RandomWalkXiFinder();
    final double[] xi = xiFinder.maximizeSingleSample(Collections.singletonList(new int[] {2, 3}), Collections.singletonList(new double[][] {{0, 0, 10, 10}}), 0).mXi;
    assertEquals(1, xi[0] + xi[1], 1e-10);
    assertEquals(0.5, xi[0], 0.1);
    assertEquals(0.5, xi[1], 0.1);
  }

  public void testEasyUnbalanced() {
    final RandomWalkXiFinder xiFinder = new RandomWalkXiFinder();
    final double[] xi = xiFinder.maximizeSingleSample(AT, Collections.singletonList(new double[][] {{10, 20, 0, 0}}), 0).mXi;
    assertEquals(1, xi[0] + xi[1], 1e-10);
    assertEquals(0.3333, xi[0], 0.1);
    assertEquals(0.6667, xi[1], 0.1);
  }

  public void test3StrainsBalanced() {
    final RandomWalkXiFinder xiFinder = new RandomWalkXiFinder();
    final double[] xi = xiFinder.maximizeSingleSample(Collections.singletonList(new int[] {0, 1, 1}), Collections.singletonList(new double[][] {{10, 10, 0, 0}}), 0).mXi;
    //System.out.println(java.util.Arrays.toString(xi));
    assertEquals(1, xi[0] + xi[1] + xi[2], 1e-10);
    assertEquals(0.5, xi[0], 0.1);
    assertEquals(0.5, xi[1] + xi[2], 0.1);
  }

  public void test4StrainsBalanced() {
    final RandomWalkXiFinder xiFinder = new RandomWalkXiFinder();
    final double[] xi = xiFinder.maximizeSingleSample(Collections.singletonList(new int[] {0, 1, 1, 0}), Collections.singletonList(new double[][] {{10, 10, 0, 0}}), 0).mXi;
    //System.out.println(java.util.Arrays.toString(xi));
    assertEquals(1, xi[0] + xi[1] + xi[2] + xi[3], 1e-10);
    assertEquals(0.5, xi[0] + xi[3], 0.1);
    assertEquals(0.5, xi[1] + xi[2], 0.1);
  }

  public void testNoisyBalanced() {
    final RandomWalkXiFinder xiFinder = new RandomWalkXiFinder();
    final double[] xi = xiFinder.maximizeSingleSample(AT, Collections.singletonList(new double[][] {{50, 50, 2, 0}}), 0).mXi;
    //System.out.println(java.util.Arrays.toString(xi));
    assertEquals(1, xi[0] + xi[1], 1e-10);
    assertEquals(0.5, xi[0], 0.1);
    assertEquals(0.5, xi[1], 0.1);
  }

  public void testThreeAlleles() {
    final RandomWalkXiFinder xiFinder = new RandomWalkXiFinder();
    final double[] xi = xiFinder.maximizeSingleSample(Collections.singletonList(new int[] {0, 1, 2}), Collections.singletonList(new double[][] {{50, 50, 50, 2}}), 0).mXi;
    //System.out.println(java.util.Arrays.toString(xi));
    assertEquals(1, xi[0] + xi[1] + xi[2], 1e-10);
    assertEquals(0.333, xi[0], 0.1);
    assertEquals(0.333, xi[1], 0.1);
    assertEquals(0.333, xi[2], 0.1);
  }

  public void testThreePositionEasy() {
    final RandomWalkXiFinder xiFinder = new RandomWalkXiFinder();
    final double[][] ev = {{5, 5, 5, 0}};
    final double[] xi = xiFinder.maximizeSingleSample(Arrays.asList(new int[] {0, 1, 2}, new int[] {2, 1, 0}, new int[] {1, 2, 0}), Arrays.asList(ev, ev, ev), 0).mXi;
    //System.out.println(java.util.Arrays.toString(xi));
    assertEquals(1, xi[0] + xi[1] + xi[2], 1e-10);
    assertEquals(0.333, xi[0], 0.1);
    assertEquals(0.333, xi[1], 0.1);
    assertEquals(0.333, xi[2], 0.1);
  }

  private void checkThreePositionHardLogArithmetic(final PossibilityArithmetic arith) {
    final RandomWalkXiFinder xiFinder = new RandomWalkXiFinder(arith);
    final double[][] ev = {{10, 10, 0, 0}};
    final List<int[]> alpha = Arrays.asList(new int[]{0, 0, 1}, new int[]{0, 1, 0}, new int[]{1, 0, 0});
    final List<double[][]> evidence = Arrays.asList(ev, ev, ev);
    final double[] xi = xiFinder.maximizeSingleSample(alpha, evidence, 0).mXi;
    //System.out.println(java.util.Arrays.toString(xi));
    final double expected = arith.prob2Poss(1.0 / 3.0);
    assertEquals(expected, xi[0], 0.1);
    assertEquals(expected, xi[1], 0.1);
    assertEquals(expected, xi[2], 0.1);
    final double[][] xiSamples = xiFinder.maximize(alpha, evidence);
    assertEquals(expected, xiSamples[0][0], 0.1);
    assertEquals(expected, xiSamples[0][1], 0.1);
    assertEquals(expected, xiSamples[0][2], 0.1);
  }
  public void test2Samples() {
    final PossibilityArithmetic arith = LogPossibility.SINGLETON;
    final RandomWalkXiFinder xiFinder = new RandomWalkXiFinder(arith);
    final List<int[]> alpha = Arrays.asList(new int[]{0, 0, 1}, new int[]{0, 1, 0}, new int[]{1, 0, 0});
    final List<double[][]> evidence = new ArrayList<>();
    evidence.add(new double[][] {{20, 10, 0, 0}, {30, 10, 0, 0}});
    evidence.add(new double[][] {{20, 10, 0, 0}, {20, 20, 0, 0}});
    evidence.add(new double[][] {{20, 10, 0, 0}, {30, 10, 0, 0}});
    final double expected = arith.prob2Poss(1.0 / 3.0);
    final double[][] xiSamples = xiFinder.maximize(alpha, evidence);
    assertEquals(expected, xiSamples[0][0], 0.1);
    assertEquals(expected, xiSamples[0][1], 0.1);
    assertEquals(expected, xiSamples[0][2], 0.1);

    assertEquals(arith.prob2Poss(0.25), xiSamples[1][0], 0.1);
    assertEquals(arith.prob2Poss(0.5), xiSamples[1][1], 0.1);
    assertEquals(arith.prob2Poss(0.25), xiSamples[1][2], 0.1);
  }

  public void testThreePositionHard() {
    checkThreePositionHardLogArithmetic(SimplePossibility.SINGLETON);
    checkThreePositionHardLogArithmetic(LogPossibility.SINGLETON);
  }

}
