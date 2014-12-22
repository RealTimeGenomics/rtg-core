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

package com.rtg.variant.util.arithmetic;

import java.io.PrintStream;

import com.rtg.util.Utils;

/**
 */
public final class Timing {

  private Timing() { }

  private static void time(final PossibilityArithmetic arith, final PrintStream ps) {
    final long iter = 100000000;
    final long t0 = System.nanoTime();
    double x = arith.one();
    final double m1 = arith.prob2Poss(0.2);
    final double m2 = arith.prob2Poss(0.99);
    final double m3 = arith.prob2Poss(0.01);
    for (long l = 0; l < iter; l++) {
      final double y1 = arith.multiply(x, m1);
      final double y2 = arith.multiply(x, m2);
      final double y3 = arith.multiply(x, m3);
      final double z1 = arith.add(y1, y2);
      x = arith.add(y3, z1);
    }
    final long t1 = System.nanoTime();
    final long delta = t1 - t0;
    final double t = delta / (double) iter;
    ps.println(arith.toString() + " " + Utils.realFormat(t, 1) + "ns");
  }

  /**
   * @param args command line arguments ignored.
   */
  public static void main(String[] args) {
    for (int i = 0; i < 3; i++) {
      //time(IntegerPossibility.SINGLETON, System.err);
      //time(SimplePossibility.SINGLETON, System.err);
      //time(LogApproximatePossibility.SINGLETON, System.err);
      time(LogPossibility.SINGLETON, System.err);
    }
  }

}
