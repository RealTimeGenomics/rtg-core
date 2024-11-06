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

package com.rtg.variant.bayes.multisample.population;


import com.rtg.util.Utils;
import com.rtg.variant.bayes.snp.DescriptionCommon;

/**
 * Plot performance of EM algorithm for YZ case where neither Y nor Z is reference and probability
 * is varied from 0 to 1.
 *
 * @param <D> description type
 */
public final class ConvergenceYZ<D extends DescriptionCommon> {

  private ConvergenceYZ() { }

  private static double itEstimators(final int samples, double errorRate, ChrType type, final Estimator estimator, final double p) {
    final ConvergencePlot convergencePlot = new ConvergencePlot(samples, type, new double[] {0.0, p, 1.0 - p, 0.0}, errorRate, estimator, null);
    return convergencePlot.expectedCoverage(20, 20000);
  }

  /**
   * @param args command line arguments. - samples, type
   */
  public static void main(String[] args) {
    final int samples = Integer.parseInt(args[0]);
    final ChrType type = ChrType.valueOf(args[1]);
    System.out.println("samples= " + samples);
    System.out.println("p   Null    HW");
    final double errorRate = 0.01;
    final int range = 10;
    for (int i = 0; i <= range; ++i) {
      final double p = (double) i / range;
      System.out.print(Utils.realFormat(p, 3));
      final double nullAvg = itEstimators(samples, errorRate, type, new NullEstimator(), p);
      System.out.print("  " + Utils.realFormat(nullAvg, 3));
      final double hwAvg = itEstimators(samples, errorRate, type, new HwEstimator(), p);
      System.out.print("  " + Utils.realFormat(hwAvg, 3));
      System.out.println();
    }
  }
}
