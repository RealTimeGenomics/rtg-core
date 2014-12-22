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


import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import com.rtg.util.Utils;
import com.rtg.variant.GenomePriorParams;
import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.Model;
import com.rtg.variant.bayes.ModelInterface;
import com.rtg.variant.bayes.multisample.population.Convergence.SimulationResult;
import com.rtg.variant.bayes.snp.HypothesesPrior;
import com.rtg.variant.bayes.snp.HypothesesSnp;
import com.rtg.variant.bayes.snp.StatisticsSnp;
import com.rtg.variant.util.arithmetic.LogPossibility;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;

/**
 * Write various outputs, suitable for plotting, of simulations of the EM algorithm.
 *
 */
public class ConvergencePlot {

  private final int mSamples;

  private final ChrType mType;

  private final double[] mDistr;

  private final double mErrorRate;

  private final Estimator mEstimator;

  private final PrintStream mOut;

  private final HypothesesPrior<Description> mDiploid;

  private final HypothesesPrior<Description> mHaploid;

  private final Random mRandom;


  /**
   * @param samples number of samples.
   * @param type of chromosome being considered.
   * @param distr probability distribution (over haploid nucleotides).
   * @param errorRate error rate for individual nucleotides.
   * @param estimator estimator to use
   * @param out output will go to this stream
   */
  ConvergencePlot(final int samples, final ChrType type, final double[] distr, final double errorRate, final Estimator estimator, final PrintStream out) {
    mSamples = samples;
    mType = type;
    mDistr = distr;
    mErrorRate = errorRate;
    mEstimator = estimator;
    mOut = out;
    //mArith = SimplePossibility.SINGLETON;
    final PossibilityArithmetic arith = LogPossibility.SINGLETON;
    final GenomePriorParams params = GenomePriorParams.builder().create();
    mDiploid = new HypothesesSnp(arith, params, false, 0);
    mHaploid = new HypothesesSnp(arith, params, true, 0);
    mRandom = new Random(143);
  }

  Convergence cleanConvergence() {
    final List<ModelInterface<?>> models = new ArrayList<>();
    for (int i = 0; i < mSamples; i++) {
      final HypothesesPrior<Description> hyp;
      if (mType == ChrType.AUTOSOMAL) {
        hyp = mDiploid;
      } else if (mType == ChrType.Y) {
        hyp = mHaploid;
      } else {
        hyp = mRandom.nextInt(2) == 0 ? mHaploid : mDiploid;
      }
      models.add(new Model<>(hyp,  new StatisticsSnp(hyp.description())));
    }
    return new Convergence(mHaploid, mDiploid, mDistr, mErrorRate, mEstimator, models, mRandom);
  }
  /**
   * Plot count of incorrect calls against a single iteration of the coverage.
   * Gives each line as:
   *  coverage
   *  total-incorrect
   *  incorrect for each nucleotide: A, C, G, T.
   * @param maxCoverage maximum coverage.
   */
  void varyCoverage(final int maxCoverage) {
    final Convergence conv = cleanConvergence();
    for (int i = 1; i <= maxCoverage; i++) {
      final SimulationResult iterate = conv.simulate();
      mOut.print(i + " " + iterate.totalIncorrect() + " ");
      final int[] incorrect = iterate.incorrect();
      for (int anIncorrect : incorrect) {
        mOut.print(anIncorrect + " ");
      }
      mOut.println();
    }
  }

  /**
   * Plot average number of incorrect calls over many iterations of the simulation.
   * Gives each line as:
   *  coverage
   *  average-incorrect
   * @param maxCoverage maximum coverage.
   * @param totalCalls totalnumber of calls counted across iterations and samples (keeps accuracy and time comparable as number of samples varied)
   * @return total ratio of incorrect calls
   */
  double iterate(final int maxCoverage, final int totalCalls) {
    final int iterations = (totalCalls + mSamples - 1) / mSamples;
    final int[] totals = new int[maxCoverage];
    double total = 0.0;
    for (int i = 0; i < iterations; i++) {
      final Convergence conv = cleanConvergence();
      for (int j = 1; j <= maxCoverage; j++) {
        final SimulationResult simres = conv.simulate();
        //System.err.println("i=" + i + " j=" + j + " incorr=" + simres.totalIncorrect());
        final int totalIncorrect = simres.totalIncorrect();
        totals[j - 1] += totalIncorrect;
        total += totalIncorrect;
      }
    }
    final double n = iterations * mSamples;
    for (int i = 0; i < maxCoverage; i++) {
      final double d = totals[i] / n;
      mOut.println((i + 1) + " " + Utils.realFormat(d, 3));
    }
    return total / (iterations * mSamples);
  }

  /**
   * Compute the average coverage at which errors occur by varying coverage and using multiple iterations.
   * @param maxCoverage maximum coverage.
   * @param iterations number of iterations.
   * @return expected coverage for an error.
   */
  double expectedCoverage(final int maxCoverage, final int iterations) {
    double total = 0.0;
    for (int i = 0; i < iterations; i++) {
      final Convergence conv = cleanConvergence();
      for (int j = 1; j <= maxCoverage; j++) {
        final SimulationResult simRes = conv.simulate();
        total += simRes.totalIncorrect();
      }
    }
    return total / (iterations * mSamples);
  }

  private static void itEstimators(final int samples, final double errorRate, final double[] prob, final Estimator estimator, final PrintStream out) {
    out.println(estimator.toString());
    final ConvergencePlot convergencePlot = new ConvergencePlot(samples, ChrType.AUTOSOMAL, prob, errorRate, estimator, out);
    final double avg = convergencePlot.iterate(20, 20000);
    final String avgStr = Utils.realFormat(avg, 3);
    out.println("Average= " + avgStr);
    out.println();
  }

  /**
   * @param args command line arguments. - samples, x, y, z (allele probability)
   */
  public static void main(String[] args) {
    final int samples = Integer.parseInt(args[0]);
    final double x = Double.parseDouble(args[2]);
    final double y = Double.parseDouble(args[2]);
    final double z = Double.parseDouble(args[2]);
    final double[] prob = {x, y, z, 0.0};
    final double errorRate = 0.05;
    System.out.println("samples= " + samples);
    itEstimators(samples, errorRate, prob, new NullEstimator(), System.out);
    itEstimators(samples, errorRate, prob, new HwEstimator(), System.out);
  }
}
