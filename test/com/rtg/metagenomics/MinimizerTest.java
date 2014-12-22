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

import junit.framework.TestCase;

/**
 */
public class MinimizerTest extends TestCase {
  Minimizer getSolver() {
    return new Minimizer();
  }

  private static class Quadratic extends Line {
    private final double mConstant;
    private final double mRoot;
    public Quadratic(final double constant, final double root) {
      mConstant = constant;
      mRoot = root;
    }

    @Override
    public double value(final double delta) {
      return values(delta)[0];
    }

    @Override
    public double[] values(double delta) {
      final double t = delta - mRoot;
      return new double[] {t * t + mConstant, 2 * t};
    }

    @Override
    public int derivativeOrder() {
      return 1;
    }
  }

  public void testQuad() {
    final Minimizer solver = getSolver();
    assertEquals(1.0, solver.solveLine(new Quadratic(1.0, 1.0), 0.0), 1e-5);
    assertEquals(0.0, solver.solveLine(new Quadratic(0.0, 0.0), 0.0), 1e-5);

    for (double k = 0; k < 10; k++) {
      assertEquals(k, solver.solveLine(new Quadratic(0.0, k), 1e-10), 1e-5);
    }
  }

  public void testInfinitePos() {
    assertEquals(0.0, getSolver().solveLine(new Quadratic(Double.POSITIVE_INFINITY, 1.0), 1e-10), 1e-5);
  }
  public void testInfiniteNeg() {
    assertEquals(0.0, getSolver().solveLine(new Quadratic(Double.NEGATIVE_INFINITY, 1.0), 1e-10), 1e-5);
  }

  public void testNaN() {
    assertEquals(0.0, getSolver().solveLine(new Quadratic(Double.NaN, 1.0), 1e-10), 1e-5);
  }


  //Similar to the case that motivated this class
  private class Complex extends Line {
    public Complex() {
    }
    @Override
    public double value(final double delta) {
      return values(delta)[0];
    }

    @Override
    public double[] values(double delta) {
      final double d2 = delta * delta;
      final double d3 = d2 * delta;
      final double d4 = d3 * delta;
      final double v0 = 0.25 * d4  - 0.5 * d3 + 0.25 * d2 + 0.01 * delta + 1.0;
      final double v1 = delta * (delta - 0.5) * (delta - 1.0) + 0.01;
      return new double[] {v0, v1};
    }

    @Override
    public int derivativeOrder() {
      return 1;
    }
  }

  public void testComplex() {
    assertEquals(0.0, getSolver().solveLine(new Complex(), 1e-10), 1e-5);
  }

}
