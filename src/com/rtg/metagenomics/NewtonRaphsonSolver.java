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
 * Root finder using Newton-Raphson method.
 */
public class NewtonRaphsonSolver extends AbstractSolver {

  /**
   * A solver with verbose output.
   * @param veryVerbose verbose output
   */
  public NewtonRaphsonSolver(boolean veryVerbose) {
    super(veryVerbose);
  }

  @Override
  protected double solveLine(Line line, double lo, double hi, double[] vlos, double[] vhis, double relThreshold) {
    double xlo = lo;
    double xhi = hi;
    //double[] xvlos = vlos;
    //double[] xvhis = vhis;

    double tlo = newtonRaphson(xlo, vlos);
    double thi = newtonRaphson(xhi, vhis);

    //System.err.println("NEW");
    while (true) {
      assert checkLineValueLo(line, xlo);
      assert checkLineValueHi(line, xlo, xhi);

      if (mVeryVerbose) {
        Diagnostic.developerLog(" lo: " + xlo + " hi: " + xhi);
      }

      if ((xhi - xlo) / (xhi + xlo) < relThreshold) {
        return (xlo + xhi) / 2.0;
      }
      //System.err.println("---------");
      //System.err.println("lo: " + xlo + " tlo: " + tlo + " vlos: " + Arrays.toString(xvlos));
      //System.err.println("hi: " + xhi + " thi: " + thi + " vhis: " + Arrays.toString(xvhis));

      double t;
      if (tlo < xhi && isFinite(tlo)) {
        if (thi > xlo && isFinite(thi)) {
          // both tlo and thi are valid
          // move both hi and lo in...
          t = (tlo + thi) / 2;
          //System.err.println("bisect tlo and thi t: " + t);
        } else {
          // thi is out of bounds
          t = tlo;
          //System.err.println("thi oob t: " + t);
        }
      } else if (isFinite(thi) && thi > xlo) {
        // tlo is out of bounds
        t = thi;
        //System.err.println("tlo oob t: " + t);
      } else {
        t = (xlo + xhi) / 2.0;  // bisect
        //System.err.println("bisect t: " + t);
      }
      if (t <= xlo || t >= xhi) {
        t = (xlo + xhi) / 2;
      }
      final double[] tvs = line.values(t);
      final double tnr = newtonRaphson(t, tvs);
      //System.err.println("t: " + t + " tvs: " + Arrays.toString(tvs) + " tnr: " + tnr);

      final double tv = tvs[0];
      if (tv > 0.0 || Double.isNaN(tv)) {
        xhi = t;
        //xvhis = tvs;
        thi = tnr;
      } else if (tv == 0.0) {
        return t;
      } else {
        xlo = t;
        //xvlos = tvs;
        tlo = tnr;
      }
    }
  }

  private static double newtonRaphson(double xn, double[] fs) {
    return xn - fs[0] / fs[1];
  }
}
