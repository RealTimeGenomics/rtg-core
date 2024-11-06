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
