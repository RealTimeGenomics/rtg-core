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
