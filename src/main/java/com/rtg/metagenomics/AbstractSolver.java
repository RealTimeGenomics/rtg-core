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
  @Override
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
    return Double.isFinite(d);
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
