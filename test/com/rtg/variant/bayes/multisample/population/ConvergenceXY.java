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

import com.rtg.util.Utils;
import com.rtg.variant.bayes.snp.DescriptionCommon;

/**
 * Plot performance of EM algorithm for XY case where X is reference and probability of X
 * is varied from 0 to 1.
 *
 * @param <D> description type
 */
public final class ConvergenceXY<D extends DescriptionCommon> {

  private ConvergenceXY() { }

  private static double itEstimators(final int samples, double errorRate, ChrType type, final Estimator estimator, final double p) {
    final ConvergencePlot convergencePlot = new ConvergencePlot(samples, type, new double[] {p, 1.0 - p, 0.0, 0.0}, errorRate, estimator, null);
    return convergencePlot.expectedCoverage(20, 20000);
  }

  /**
   * @param args command line arguments. - samples
   */
  public static void main(String[] args) {
    final int sampless = Integer.parseInt(args[0]);
    final ChrType type = ChrType.valueOf(args[1]);
    System.out.println("samples= " + sampless);
    System.out.println("p   Null    HW");
    final double errorRate = 0.01;
    final int range = 10;
    for (int i = 0; i <= range; i++) {
      final double p = ((double) i) / range;
      System.out.print(Utils.realFormat(p, 3));
      final double nullAvg = itEstimators(sampless, errorRate, type, new NullEstimator(), p);
      System.out.print("  " + Utils.realFormat(nullAvg, 3));
      final double hwAvg = itEstimators(sampless, errorRate, type, new HwEstimator(), p);
      System.out.print("  " + Utils.realFormat(hwAvg, 3));
      System.out.println();
    }
  }
}
