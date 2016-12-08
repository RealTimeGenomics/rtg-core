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

import java.util.Arrays;

import com.rtg.variant.util.arithmetic.LogApproximatePossibility;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;

import junit.framework.TestCase;

/**
 */
public class AlphaSelectorTest extends TestCase {

  public void test() {
    final double[][] xi = {{0.8, 0.2}, {0.2, 0.8}};
    final double[][] reads = {{2, 8, 0, 0}, {8, 2, 0, 0}};
    final int[] expected = {1, 0};
    check(xi, reads, expected);
  }
  public void testInexact() {
    final double[][] xi = {{0.8, 0.2}, {0.2, 0.8}};
    final double[][] reads = {{0, 0, 1, 9}, {0, 0, 7, 3}};
    final int[] expected = {3, 2};
    check(xi, reads, expected);
  }

  public void testNotFirst() {
    final double[][] xi = {{0.8, 0.2}, {0.2, 0.8}};
    final double[][] reads = {{0, 0, 1, 9}, {0, 0, 7, 3}};
    final int[] expected = {3, 2};
    check(xi, reads, expected);
  }

  public void test3strains() {
    final double[][] xi = {{0.5, 0.2, 0.3}, {0.1, 0.8, 0.1}};
    // for reads we would use {{0, 0, 2, 8}, {0, 0, 8, 2}} except there isn't quite enough evidence to override the beta.
    //
    final double[][] reads = {{0, 0, 2, 10}, {0, 0, 8, 4}};
    final int[] expected = {3, 2, 3};
    check(xi, reads, expected);
  }

  public void test3samples() {
    final double[][] xi = {{0.5, 0.2, 0.3}, {0.1, 0.8, 0.1}, {0.1, 0.3, 0.6}};
    final double[][] reads = {{0, 0, 3, 7}, {0, 0, 9, 1}, {0, 0, 4, 6}};
    final int[] expected = {3, 2, 3};
    check(xi, reads, expected);
  }
  public void testRefIsLikely() {
    final double[][] xi = {{0.5, 0.2, 0.3}, {0.1, 0.8, 0.1}, {0.1, 0.3, 0.6}};
    final double[][] reads = {{0, 0, 0, 1}, {0, 0, 1, 0}, {2, 0, 0, 0}};
    final int[] expected = {0, 0, 0};
    check(xi, reads, expected);

  }

  private void check(double[][] xi, double[][] reads, int[] expected) {
    final AlphaScore alphaScore = getAlphaScore(xi, reads);
    final int[] assignments = alphaScore.mCalls;
    assertTrue("expected <" + Arrays.toString(expected) + "> but was <" + Arrays.toString(assignments) + ">", Arrays.equals(expected, assignments));
  }

  private AlphaScore getAlphaScore(double[][] xi, double[][] reads) {
    final int template = 0;
    final double[] beta = new double[xi[0].length];
    Arrays.fill(beta, 0.001);
    final double error = 0.01;
    final PossibilityArithmetic arith = LogApproximatePossibility.SINGLETON;
    final double[][] thetaLookup = AlphaSelector.computeThetaLookup(xi, arith, arith.prob2Poss(1 - error), arith.prob2Poss(error / 3));
    return AlphaSelector.alphaPosition(template, beta, reads, thetaLookup, arith, xi[0].length);
  }

  public void testScoreConvergence() {
    // Test that evidence closer to the expected ratio gives higher scores
    final double[][] xi = {{0.5, 0.2, 0.3}, {0.1, 0.8, 0.1}, {0.1, 0.3, 0.6}};
    final double[][] reads = {{0, 0, 3, 7}, {0, 0, 9, 1}, {0, 0, 4, 6}};
    final double[][][] betterReads = {
        {{0, 0, 3, 8}, {0, 0, 9, 1}, {0, 0, 4, 6}}
        , {{0, 0, 3, 7}, {0, 0, 9, 2}, {0, 0, 4, 6}}
        , {{0, 0, 3, 7}, {0, 0, 9, 1}, {0, 0, 4, 7}}
    };
    final double initialScore = getAlphaScore(xi, reads).mLikelihood;
    for (final double[][] betterRead : betterReads) {
      final double score = getAlphaScore(xi, betterRead).mLikelihood;
      assertTrue("Not better: " + Arrays.deepToString(betterRead), score > initialScore);
    }
    final double[][][] worseReads = {
        {{0, 0, 4, 7}, {0, 0, 9, 1}, {0, 0, 4, 6}}
        , {{0, 0, 3, 7}, {0, 0, 10, 1}, {0, 0, 4, 6}}
        , {{0, 0, 3, 7}, {0, 0, 9, 1}, {0, 0, 5, 6}}
    };
    for (int i = 0; i < betterReads.length; ++i) {
      final double score = getAlphaScore(xi, worseReads[i]).mLikelihood;
      assertTrue("Not worse: " + Arrays.deepToString(worseReads[i]), score < initialScore);
    }

  }

}
