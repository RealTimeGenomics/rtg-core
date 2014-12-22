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

import com.rtg.util.Utils;
import com.rtg.variant.bayes.snp.DescriptionCommon;

/**
 * Plot performance of EM algorithm for XY case where X is reference and probability of X
 * is varied from 0 to 1.
 *
 * @param <D> description type
 */
public final class ConvergenceXYZ<D extends DescriptionCommon> {

  private ConvergenceXYZ() { }

  private static void itEstimators(final int samples, double errorRate, ChrType type, final Estimator estimator, final PrintStream out) {
    out.println(estimator.toString());
    final double third = 1.0 / 3.0;
    final ConvergencePlot convergencePlot = new ConvergencePlot(samples, type, new double[] {third, third, third, 0.0}, errorRate, estimator, out);
    final double avg = convergencePlot.iterate(20, 20000);
    final String avgStr = Utils.realFormat(avg, 3);
    out.println("Average= " + avgStr);
    out.println();
  }

  /**
   * @param args command line arguments. - samples
   */
  public static void main(String[] args) {
    final int samples = Integer.parseInt(args[0]);
    final ChrType type = ChrType.valueOf(args[1]);
    final double errorRate = 0.01;
    System.out.println("samples= " + samples);
    itEstimators(samples, errorRate, type, new NullEstimator(), System.out);
    itEstimators(samples, errorRate, type, new HwEstimator(), System.out);
  }
}
