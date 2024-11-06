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
    Quadratic(final double constant, final double root) {
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

    for (double k = 0; k < 10; ++k) {
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
    Complex() {
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
