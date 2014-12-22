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
