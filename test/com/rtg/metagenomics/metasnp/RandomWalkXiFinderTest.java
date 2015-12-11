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
