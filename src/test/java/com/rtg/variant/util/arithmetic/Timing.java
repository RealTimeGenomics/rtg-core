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
    for (long l = 0; l < iter; ++l) {
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
    for (int i = 0; i < 3; ++i) {
      //time(IntegerPossibility.SINGLETON, System.err);
      //time(SimplePossibility.SINGLETON, System.err);
      //time(LogApproximatePossibility.SINGLETON, System.err);
      time(LogPossibility.SINGLETON, System.err);
    }
  }

}
