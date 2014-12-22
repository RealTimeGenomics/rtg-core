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
import java.util.Iterator;
import java.util.List;
import java.util.Random;

import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;
import com.rtg.variant.util.arithmetic.SimplePossibility;

/**
 * Monte Carlo solver.  Actually closer to a hill climb than Monte Carlo.
 * See "Strain Detection in Metagenomic Samples" by John G. Cleary.
 */
public class RandomWalkXiFinder {
  
  private final Random mRandom = new Random(42);
  private static final int ITERATIONS = 1000; // What is a sensible number for this?
  private static final int MONTE_CARLO_ITERATIONS = 10; // 1 == do a single random start point
  
  private final PossibilityArithmetic mArith;
  
  RandomWalkXiFinder(final PossibilityArithmetic arith) {
    mArith = arith;
  }
  
  RandomWalkXiFinder() {
    this(SimplePossibility.SINGLETON);
  }
  
  static void normalize(final double[] p, final PossibilityArithmetic arith) {
    // In situ arithmetic aware normalize.  Assumes entries are positive.
    double t = arith.zero();
    for (final double v : p) {
      t = arith.add(t, v);
    }
    for (int k = 0; k < p.length; k++) {
      p[k] = arith.divide(p[k], t);
    }
  }
  
  private double[] perturb(final double[] xi) {
    final int entry = mRandom.nextInt(xi.length);
    final double[] modifiedXi = Arrays.copyOf(xi, xi.length);
    final double perturbation = mArith.prob2Poss(mRandom.nextDouble());
    final double current = modifiedXi[entry];
    if (mRandom.nextBoolean()) {
      modifiedXi[entry] = mArith.divide(current, perturbation);  // increasing
    } else {
      modifiedXi[entry] = mArith.multiply(current, perturbation);  // decreasing
    }
    normalize(modifiedXi, mArith);
    return modifiedXi;
  }

  private double[] kurtTransform(final int numStrains, final List<int[]> alpha, final List<int[][]> count, int sample) {
    // We want to perform a massive product over all the positions.  Since the
    // product is commutative we can accumulate a count for each possible
    // combination of xi values.  This must be done in a way consistent with
    // the allSums procedure.
    final double[] powers = new double[1 << numStrains];
    final Iterator<int[]> alphaIterator = alpha.iterator();
    for (final int[][] cnt : count) {
      final int[] a = alphaIterator.next();
      for (int allele = 0; allele < cnt[sample].length; allele++) {
        int index = 0;
        for (int k = numStrains - 1; k >= 0; k--) {
          index <<= 1;
          if (a[k] == allele) {
            index++;
          }
        }
        powers[index] += cnt[sample][allele];
      }
    }
    return powers;
  }
  
  private double[] allSums(final double[] xi) {
    // Given an array xi of length n compute an array gamma of length 2^n
    // comprising all possible sums of the elements of xi.  In particular, if
    // k = sum_i b_i2^i is the binary representation of k, then
    // gamma[k] = sum_{b[i]=1} xi[i].  Note gamma[0] is always 0 and
    // gamma[2^n-1] is the sum of all the elements of xi.  Also,
    // gamma[2^i]=xi[i].  Uses the arithmetic of the class.
    final double[] sums = new double[1 << xi.length];
    Arrays.fill(sums, mArith.zero());
    int j = 1;
    for (final double x : xi) {
      for (int i = 0; i < j; i++) {
        sums[j + i] = mArith.add(sums[i], x);
      }
      j <<= 1;
    }
    return sums;
  }

  static double[] uniformDistribution(final int numStrains, PossibilityArithmetic arith) {
    // Construct a distribution where each strain is equally likely
    // Returned array is in the arithmetic of the class
    final double[] xi = new double[numStrains];
    Arrays.fill(xi, arith.prob2Poss(1.0 / numStrains));
    return xi;
  }

  XiScore maximizeSingleSample(final List<int[]> alpha, final List<int[][]> evidence, int sample) {
    // Let Xi=(xi_0,...,xi_n) satisfy 0 <= xi[i] <= 1, sum x[i]=1.
    // Maximize product_x product_a theta_x[a]^count_x[a] where
    //   theta_x[a] = sum_{alpha_x(i)=a} xi[i]

    assert alpha.size() == evidence.size(); // indexed over positions x
    
    final int numStrains = alpha.iterator().next().length; // Ugly
    final double[] powers = kurtTransform(numStrains, alpha, evidence, sample);

    double best = mArith.zero();
    double[] bestXi = null;
    
    for (int monteCarloIteration = 0; monteCarloIteration < MONTE_CARLO_ITERATIONS; monteCarloIteration++) {
      // Start from a uniform assignment.  It should be possible to get a
      // better initial estimate.  There did not appear to be any advantage
      // starting from a random distribution.
      double[] xi = uniformDistribution(numStrains, mArith);

      // Hill climbing by random perturbation of current solution.
      for (int hillClimbIteration = 0; hillClimbIteration < ITERATIONS; hillClimbIteration++) {
        final double[] newXi = perturb(xi);
        final double[] xiSums = allSums(newXi);
        double p = mArith.one();
        for (int j = 1; j < xiSums.length; j++) {
          if (mArith.gt(xiSums[j],  mArith.zero())) {
            final double r = mArith.pow(xiSums[j], powers[j]);
            p = mArith.multiply(p, r);
          }
        }
        if (mArith.gt(p, best)) {
          best = p;
          xi = newXi;
          bestXi = xi;
        }
      }
    }
    Diagnostic.userLog("lnP(A|alpha,Xi): " + mArith.poss2Ln(best) + " sample: " + sample);
    return new XiScore(best, bestXi);
  }

  /**
   * Compute a best estimate for the ratio of species in each sample.
   * The first parameter, alpha, has for each position an assignment of alleles to
   * strains. That it, each of the integer arrays in alpha has length equal to the
   * number of strains and contains the putative allele for that strain for a position.
   * The outer list of the evidence parameter is over position. The outer array is
   * over sample.  Each of the individual integer arrays give the allele counts
   * observed in the reads for the given position and sample.  Typically this will
   * be four nucleotides long, representing the number of A, C, G, and T nucleotides
   * in the reads.
   * @param alpha allele assignments
   * @param evidence allele counts as evidenced in the reads
   * @return best estimate for xi in arithmetic
   */
  double[][] maximize(final List<int[]> alpha, final List<int[][]> evidence) {
    // John's theory says we can maximize each sample independently during this step.
    final double[][] res = new double[evidence.get(0).length][];
    double p = mArith.one();
    for (int sample = 0; sample < evidence.get(0).length; sample++) {
      final XiScore xi = maximizeSingleSample(alpha, evidence, sample);
      res[sample] = xi.mXi;
      p = mArith.multiply(p, xi.mScore);
    }
    Diagnostic.userLog("-lnP(Xi|alpha,A): " + -mArith.poss2Ln(p));
    return res;
  }

}
