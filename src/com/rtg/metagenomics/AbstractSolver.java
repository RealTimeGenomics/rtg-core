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
 * Looks for a zero in the function.
 */
public abstract class AbstractSolver implements LineSolver {
  final boolean mVeryVerbose;

  /**
   * A solver with verbose output.
   * @param veryVerbose verbose output
   */
  AbstractSolver(boolean veryVerbose) {
    mVeryVerbose = veryVerbose;
  }

  /**
   * Finds the point between <code>lo</code> and <code>hi</code> where the <code>line</code> crossed 0.
   * It is expected that the line crosses 0 within the given range.
   * @param line line the process
   * @param lo low x bound of region to search
   * @param hi high x bound of region to search
   * @param vlos optional array of y values and derivatives for the low end of the range
   * @param vhis optional array of y values and derivatives for the high end of the range
   * @param relThreshold termination threshold
   * @return cross over point
   */
  protected abstract double solveLine(Line line, double lo, double hi, double[] vlos, double[] vhis, double relThreshold);

  /**
   * Finds the point between <code>lo</code> and <code>hi</code> where the <code>line</code> crossed 0.
   * It is expected that the line crosses 0 within the given range.
   * @param line line the process
   * @param lo low x bound of region to search
   * @param hi high x bound of region to search
   * @param relThreshold termination threshold
   * @return cross over point
   */
  protected double solveLine(Line line, double lo, double hi, double relThreshold) {
    return solveLine(line, lo, hi, line.values(lo), line.values(hi), relThreshold);
  }

  /**
   * Finds point where line crosses y at 0.  Assumes this is at a point where x is greater than 0.
   * @param line line the process
   * @param relThreshold termination threshold
   * @return cross over point
   */
  public double solveLine(Line line, double relThreshold) {
    if (mVeryVerbose) {
      Diagnostic.developerLog("starting solveLine");
    }

    double lo = 0.0;
    double[] vlos = line.values(lo);
    final double v0 = vlos[0];

    if (v0 >= 0.0 || Double.isNaN(v0)) {
      Diagnostic.developerLog("Initial value v0: " + v0);
      return 0.0;
    }
    assert v0 < 0.0 && !Double.isNaN(v0) : v0;

    double hi = 1.0;
    double[] vhis;
    //keep doubling till infinite or evaluates to positive
    while (true) {
      assert isFinite(lo) && 0.0 <= lo;
      assert isFinite(hi) && lo < hi;
      assert checkLineValueLo(line, lo);
      if (mVeryVerbose) {
        Diagnostic.developerLog(" lo: " + lo + " hi: " + hi);
      }
      vhis = line.values(hi);
      final double vhi = vhis[0];
      if (Double.isNaN(vhi)) {
        if (mVeryVerbose) {
          Diagnostic.developerLog(" NaN when evaluating hi: " + hi);
        }
        break; //but in the spirit of robustness we continue
      }
      if (vhi > 0.0) {
        break;
      }
      if (vhi == 0.0) {
        return hi;
      }
      final double h2 = hi * 2.0;
      if (Double.isInfinite(h2)) {
        Diagnostic.developerLog("Infinite case. Defaulting to 1.0");
        return 1.0; //there is some controversy about this
      }
      lo = hi;
      vlos = vhis;
      hi = h2;
    }
    return solveLine(line, lo, hi, vlos, vhis, relThreshold);
  }

  /**
   * Returns whether the give number is finite.
   * @param d number to test
   * @return true if finite
   */
  static boolean isFinite(final double d) {
    return !Double.isInfinite(d) && !Double.isNaN(d);
  }

  protected static boolean checkLineValueLo(final Line line, final double lo) {
    assert lo >= 0.0 && isFinite(lo);
    final double lov = line.value(lo);
    assert !Double.isNaN(lov) && lov <= 0.0;
    return true;
  }

  protected static boolean checkLineValueHi(final Line line, final double lo, final double hi) {
    assert hi > 0.0 && hi > lo && isFinite(hi);
    final double hiv = line.value(hi);
    assert Double.isNaN(hiv) || hiv > 0.0 : hiv;
    return true;
  }
}
