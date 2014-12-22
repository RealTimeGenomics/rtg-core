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
package com.rtg.metagenomics;

import com.rtg.util.diagnostic.Diagnostic;

/**
 * Base class for different types of root finding mechanisms.
 * Looks for a minimum in the function which must be less than the values at the limits.
 */
class Minimizer {

  Minimizer() {
  }

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
