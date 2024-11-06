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
import java.util.List;

import com.rtg.variant.util.arithmetic.PossibilityArithmetic;

/**
 */
public final class AlphaSelector {

  /** Highest allele value */
  static final int MAX_VALUE = 3;

  private AlphaSelector() { }

  static void updateThetaMask(int[] masks, int strain, int base) {
    for (int i = 0; i < masks.length; ++i) {
      masks[i] &= ~(1 << strain);
    }
    masks[base] |= 1 << strain;
  }

  /**
   * Faster implementation of alpha position using a precomputed theta table and match patterns
   *
   * @param referenceAllele base in the reference at this position
   * @param probSpaceBeta probability of a variant in each strain (in prob space)
   * @param reads evidence for this position as count of reads with each base. first index is sample, second index is allele
   * @param thetaLookup precomputed theta table
   * @param arith arithmetic object
   * @param nStrains number of strains
   * @return assignments and a score.
   */
  static AlphaScore alphaPosition(int referenceAllele, double[] probSpaceBeta, double[][] reads, double[][] thetaLookup, PossibilityArithmetic arith, int nStrains) {
    return alphaPosition(referenceAllele, new ProbAlphaSimpleBeta(probSpaceBeta), reads, thetaLookup, arith, nStrains);
  }
  static AlphaScore alphaPosition(int referenceAllele, ProbAlpha pAlpha, double[][] reads, double[][] thetaLookup, PossibilityArithmetic arith, int nStrains) {
    final int[] strainVariants = new int[nStrains];
    int stackPos = 0;
    double bestScore = arith.zero();
    double bestEvidenceScore = arith.zero();
    double restScore = arith.zero();
    List<Integer> best = new ArrayList<>(nStrains);
    for (int i = 0; i < nStrains ; ++i) {
      best.add(0);
    }
    final int[] thetaMask = new int[MAX_VALUE + 1];
    int last = -1;
    // Loop over all possible alpha_x assigments (assigments of alleles to strains)
    while (stackPos > 0 || last != MAX_VALUE) {
      strainVariants[stackPos++] = last + 1;
      updateThetaMask(thetaMask, stackPos - 1, strainVariants[stackPos - 1]);

      while (stackPos < nStrains) {
        strainVariants[stackPos++] = 0;
        updateThetaMask(thetaMask, stackPos - 1, 0);
      }

      double evidenceScore = arith.one();
      for (int sampleIndex = 0; sampleIndex < reads.length; ++sampleIndex) {
        final double[] sample = reads[sampleIndex];
        for (int i = 0; i < sample.length; ++i) {
          evidenceScore = arith.multiply(evidenceScore, arith.pow(thetaLookup[sampleIndex][thetaMask[i]], sample[i]));
        }
      }
      final double alphaScore = arith.prob2Poss(pAlpha.pAlpha(referenceAllele, strainVariants));
      final double currentScore = arith.multiply(alphaScore, evidenceScore);
//      System.err.println("score=" + currentScore + " best= " + bestScore + " strainVariants=" + strainVariants);
      if (currentScore > bestScore) {
        restScore = arith.add(restScore, bestScore);
        bestScore = currentScore;
        bestEvidenceScore = evidenceScore;
        best = new ArrayList<>(strainVariants.length);
        for (int base : strainVariants) {
          best.add(base);
        }
      } else {
        restScore = arith.add(restScore, currentScore);
      }

      while (stackPos > 0 && (last = strainVariants[--stackPos]) == MAX_VALUE) {
        // nop
      }
    }
    final int[] res = new int[best.size()];
    for (int i = 0; i < res.length; ++i) {
      res[i] = best.get(i);
    }
    return new AlphaScore(arith.divide(bestScore, restScore), bestEvidenceScore, res);
  }

  /**
   * Compute an array where the first index is sample and the second index is a bit mask for the matching strains and the value is the corresponding theta
   * @param xi current xi estimate
   * @param arith arithmetic used for scoring
   * @param notError probability of no error in poss space
   * @param thirdError one third the probability of error in poss space
   * @return an theta lookup table per sample
   */
  static double[][] computeThetaLookup(double[][] xi, PossibilityArithmetic arith, double notError, double thirdError) {
    final int combinations = 1 << xi[0].length;
    final double[][] thetaLookup = new double[xi.length][combinations];
    for (int sample = 0; sample < xi.length; ++sample) {
      final double[] x = xi[sample];
      final double[] theta = thetaLookup[sample];
      Arrays.fill(theta, arith.zero());
      for (int comb = 0; comb < combinations; ++comb) {
        double oneMinus = arith.zero();
        for (int strain = 0; strain < x.length; ++strain) {
          if ((comb & (1 << strain)) != 0) {
            theta[comb] = arith.add(theta[comb], x[strain]);
          } else {
            oneMinus = arith.add(oneMinus, x[strain]);
          }
        }
        theta[comb] = arith.add(arith.multiply(notError, theta[comb]), arith.multiply(thirdError, oneMinus));
      }
    }

    return thetaLookup;
  }
}
