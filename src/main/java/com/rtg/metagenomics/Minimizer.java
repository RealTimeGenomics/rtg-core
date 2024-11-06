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
 * Base class for different types of root finding mechanisms.
 * Looks for a minimum in the function which must be less than the values at the limits.
 */
class Minimizer {

  /**
   * Finds point where line is a minimum.  Assumes this is at a point where x is greater than 0.
   * @param line the function being minimized.
   * @param relThreshold termination threshold
   * @return the minimum point.
   */
  double solveLine(Line line, double relThreshold) {
    assert line.derivativeOrder() >= 1;

    double lo = 0.0;
    double vlo = line.value(lo);

    if (Double.isNaN(vlo)) {
      Diagnostic.developerLog("Initial value vlo: " + vlo);
      return 0.0;
    }

    double hi = 1.0;
    double vhi = line.value(hi);
    //keep doubling till NaN or greater than initial lo value
    while (true) {
      //Diagnostic.developerLog("A lo=" + lo + " vlo=" + vlo + " hi=" + hi + " vhi=" + vhi);
      assert AbstractSolver.isFinite(lo) && 0.0 <= lo;
      assert AbstractSolver.isFinite(hi) && lo < hi;
      //assert checkLineValueLo(line, lo);
      if (vhi >= vlo || Double.isNaN(vhi)) {
        break;
      }

      hi = hi * 2.0;
      vhi = line.value(hi);
    }
    //keep bisecting till get to minimum.
    //TODO optimize
    while (true) {
      //Diagnostic.developerLog("B lo=" + lo + " vlo=" + vlo + " hi=" + hi + " vhi=" + vhi);
      final double t = (lo + hi) / 2.0;
      if (t <= lo || t >= hi) {
        //Diagnostic.developerLog("t=" + t);
        return t;
      }
      final double[] tvs = line.values(t);
      final double tv = tvs[0];
      if (tv >= vlo || Double.isNaN(tv)) {
        hi = t;
        vhi = tv;
        continue;
      }
      final double tv1 = tvs[1];
      //Diagnostic.developerLog("tv1=" + tv1);
      if (tv1 <= 0.0) {
        lo = t;
        vlo = tv;
      } else {
        hi = t;
        vhi = tv;
      }
    }
  }
}
