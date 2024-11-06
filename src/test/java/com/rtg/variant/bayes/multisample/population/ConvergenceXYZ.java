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
