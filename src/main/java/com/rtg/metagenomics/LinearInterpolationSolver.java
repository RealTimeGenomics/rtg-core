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
 * Finds root via repetitive bisection using linear interpolation to estimate mid point.
 */
public class LinearInterpolationSolver extends AbstractSolver {

  /**
   * Default constructor.
   */
  public LinearInterpolationSolver() {
    this(false);
  }

  /**
   * A solver with verbose output.
   * @param veryVerbose verbose output
   */
  public LinearInterpolationSolver(boolean veryVerbose) {
    super(veryVerbose);
  }

  @Override
  protected double solveLine(Line line, double lo, double hi, double[] vlos, double[] vhis, double relThreshold) {
    // bisection
    // linear interpolation of mid point
    assert vlos != null && vlos.length >= 1;
    assert vhis != null && vhis.length >= 1;

    double xlo = lo;
    double xhi = hi;
    double vlo = vlos[0];
    double vhi = vhis[0];

    boolean linear = true;

    //System.err.println("NEW lo " + xlo + " hi " + xhi + " vlo " + vlo + " vhi " + vhi);
    //int count = 0;
    while (true) {
      assert checkLineValueLo(line, xlo);
      assert checkLineValueHi(line, xlo, xhi);

      if (mVeryVerbose) {
        Diagnostic.developerLog(" lo: " + xlo + " hi: " + xhi);
      }
      //System.err.println("x lo " + xlo + " hi " + xhi + " vlo " + vlo + " vhi " + vhi);
      //      final double maxDeriv = Math.max(-vlo, vhi);
      //      if ((xhi - xlo) * maxDeriv < relThreshold) {
      //        return (xlo + xhi) / 2.0;
      //      }

      double t;
      // linear interpolation or bisection every second iteration
      // linear interpolation only if end point values are usable
      if (linear && isFinite(vlo) && isFinite(vhi)) {
        t = xlo + (xhi - xlo) * (0 - vlo) / (vhi - vlo);
      } else {
        t = (xlo + xhi) / 2.0;
      }
      if (!isFinite(t)) { // handle case where hi and lo are very close...
        t = (xlo + xhi) / 2.0;
      }

      //System.err.println("t " + t);
      if (t <= xlo || t >= xhi) {
        return t;
      }
      final double tv = line.value(t);
      //System.err.println("tv " + tv);
      if (tv > 0.0 || Double.isNaN(tv)) {
        xhi = t;
        vhi = tv;
      } else if (tv == 0.0) {
        return t;
      } else {
        xlo = t;
        vlo = tv;
      }
      linear = !linear;
    }
  }

}
