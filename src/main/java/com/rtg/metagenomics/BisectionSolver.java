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

package com.rtg.metagenomics;

import com.rtg.util.diagnostic.Diagnostic;

/**
 * Finds root by simple repetitive bisection of bounds.
 */
public class BisectionSolver extends AbstractSolver {

  /**
   * A solver with verbose output.
   * @param veryVerbose verbose output
   */
  public BisectionSolver(boolean veryVerbose) {
    super(veryVerbose);
  }

  @Override
  protected double solveLine(Line line, double lo, double hi, double[] vlos, double[] vhis, double relThreshold) {
    double xlo = lo;
    double xhi = hi;

    while (true) {
      assert checkLineValueLo(line, xlo);
      assert checkLineValueHi(line, xlo, xhi);

      if (mVeryVerbose) {
        Diagnostic.developerLog(" lo: " + xlo + " hi: " + xhi);
      }

      if ((xhi - xlo) / (xhi + xlo) < relThreshold) {
        return (xlo + xhi) / 2.0;
      }

      final double t = (xlo + xhi) / 2.0;
      if (t <= xlo || t >= xhi) {
        return t;
      }
      final double tv = line.value(t);
      if (tv > 0.0 || Double.isNaN(tv)) {
        xhi = t;
      } else if (tv == 0.0) {
        return t;
      } else {
        xlo = t;
      }
    }
  }

}
