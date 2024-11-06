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

/**
 */
public final class EmIterate {
  /**
   * Possible beta implementations
   */
  public enum BetaType {
    /** Use a static value */
    STATIC,
    /** Re-estimate variation each iteration */
    REESTIMATE,
    /** More sophisticated re estimation which allows for between species correlation */
    COMPLEX
  }
  private EmIterate() { }
  static final class EmResult {
    final double[][] mXi;
    final List<AlphaScore> mAssignments;

    EmResult(double[][] xi, List<AlphaScore> assignments) {
      mXi = xi;
      mAssignments = assignments;
    }
  }

  interface Termination {
    /**
     * @param iteration counter indicating how many iterations so far
     * @param results collection of results from previous iterations
     * @return true if iteration should finish
     */
    boolean finished(int iteration, List<EmResult> results);
  }

  /**
   * Run for a fixed number of iterations
   */
  static class FixedIterations implements Termination {
    final int mIterations;

    FixedIterations(int iterations) {
      mIterations = iterations;
    }

    @Override
    public boolean finished(int iteration, List<EmResult> results) {
      return iteration >= mIterations;
    }
  }
  private static ProbAlpha getProbAlpha(BetaType type, List<Integer> ref, final List<int[]> assignments, long length, double[] staticBeta) {
    switch (type) {
      case COMPLEX:
        return new ComplicatedBeta(ref, assignments, length);
      case REESTIMATE:
        return new ProbAlphaSimpleBeta(estimateBeta(ref, assignments, length));
      default:
        return new ProbAlphaSimpleBeta(staticBeta);
    }
  }
  
  private static double[] estimateBeta(List<Integer> ref, final List<int[]> assignments, long length) {
    assert ref.size() == assignments.size();
    assert length >= ref.size() && length > 0;
    final double[] beta = new double[assignments.get(0).length];
    for (int pos = 0; pos < ref.size(); ++pos) {
      final int refNt = ref.get(pos);
      final int[] alpha = assignments.get(pos);
      for (int strain = 0; strain < beta.length; ++strain) {
        if (alpha[strain] != refNt) {
          beta[strain]++;
        }
      }
    }
    for (int strain = 0; strain < beta.length; ++strain) {
      beta[strain] /= length;
    }
    return beta;
  }
  
  private static double codeLength(final List<AlphaScore> scores, final int strains, final PossibilityArithmetic arith, long length) {
    // Length not currently used 
    double codeLength = 0;
    for (final AlphaScore score : scores) {
      codeLength += arith.poss2Ln(score.mBestPoss);
    }
    // Cost of xi is negligible?
    // Cost of specifying A + cost of specifying alpha
    return -codeLength + strains * scores.size() * Math.log(4);
  }

  static List<EmResult> iterate(List<Integer> ref, List<double[][]> evidence, int strains, BetaType betaType, double error, double[][] xi) {
    return iterate(ref, evidence, strains, ref.size(), LogPossibility.SINGLETON, new FixedIterations(10), betaType, error, xi);
  }
  /**
   *
   * @param ref per interesting position reference byte
   * @param evidence per position evidence position/sample/base = count
   * @param strains  number of strains to estimate
   * @param arith arithmetic space to use
   * @param approxLength approximate genome length
   * @param terminate when should we stop
   * @param updateBeta should beta be recomputed each iteration
   * @param error error rate
   * @param xiPrior priors for xi
   * @return assignments and predicted xi
   */
  static List<EmResult> iterate(List<Integer> ref, List<double[][]> evidence, int strains, long approxLength, PossibilityArithmetic arith, Termination terminate, BetaType updateBeta, double error, double[][] xiPrior) {
    final double[] beta = new double[strains];
    Arrays.fill(beta, 0.001);
    double[][] xi = xiPrior;
    List<AlphaScore> scores;
    List<int[]> assignments;
    final RandomWalkXiFinder randomWalkXiFinder = new RandomWalkXiFinder(arith);
    final List<EmResult> results = new ArrayList<>();
    ProbAlpha pAlpha = getProbAlpha(BetaType.STATIC, ref, Collections.emptyList(), 0, beta);
    for (int emIterations = 0; !terminate.finished(emIterations, results); ++emIterations) {
      Diagnostic.progress("Starting Iteration: " + emIterations);
      final int n = evidence.size();
      assignments = new ArrayList<>(n);
      scores = new ArrayList<>(n);
      int percent = 0;
      int lastReport = 0;
      final double[][] outerTheta = AlphaSelector.computeThetaLookup(xi, arith, arith.prob2Poss(1 - error), arith.prob2Poss(error / 3));
      for (int x = 0; x < n; ++x) {

        final AlphaScore alphaScore =  AlphaSelector.alphaPosition(ref.get(x), pAlpha, evidence.get(x), outerTheta, arith, xi[0].length);
        scores.add(alphaScore);
        assignments.add(alphaScore.mCalls);
        if ((x * 100 / n) > percent && x > lastReport + 10000) {
          percent = x * 100 / n;
          lastReport = x;
          Diagnostic.progress("Strain assignment: " + percent + "%");
        }
      }
      Diagnostic.userLog("Alpha code length: " + codeLength(scores, strains, arith, approxLength));
      pAlpha = getProbAlpha(updateBeta, ref, assignments, approxLength, beta);
      final double[][] newXi = randomWalkXiFinder.maximize(assignments, evidence);
      Diagnostic.userLog("Frobenius: " + Frobenius.frobeniusDistance(arith, xi, newXi));
      xi = newXi;
      Diagnostic.progress("Finished Iteration: " + emIterations);
      Diagnostic.info("Iteration " + emIterations + " beta: " + pAlpha);
      results.add(new EmResult(xi, scores));
    }
    return results;
  }

}
