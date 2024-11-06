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

import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.MemoryPrintStream;

import junit.framework.TestCase;

/**
 */
public class LoggingLineSolverTest extends TestCase {


  @Override
  protected void setUp() {
    Diagnostic.setLogStream();
  }

  private LoggingLineSolver getLineSolver(boolean verbose) {
    return new LoggingLineSolver(new BisectionSolver(verbose), verbose);
  }

  private static class LinearLine extends Line {
    private final double mConstant;
    private final double mMultiplier;
    LinearLine(final double constant, final double multiplier) {
      mConstant = constant;
      mMultiplier = multiplier;
    }
    @Override
    public double value(final double delta) {
      return mMultiplier * (mConstant + delta);
    }
  }

  private static class ZeroLine extends Line {
    @Override
    public double value(final double delta) {
      if (delta < 2.7) {
        return -1.0;
      }
      return 0.0;
    }
  }


  public void testLinear() {
    for (double k = 0; k < 10; ++k) {
      assertEquals(k, getLineSolver(false).solveLine(new LinearLine(-k, 1.0), 1e-10));
    }
  }

  public void testLinearSmall() {
    for (double k = 0; k < 10; ++k) {
      assertEquals(k, getLineSolver(false).solveLine(new LinearLine(-k, 0.5 / (k + 1)), 1e-10));
    }
  }

  public void testInfinitePos() {
    assertEquals(0.0, getLineSolver(false).solveLine(new LinearLine(Double.POSITIVE_INFINITY, 1.0), 1e-10));
  }

  public void testInfiniteNeg() {
    final MemoryPrintStream ps = new MemoryPrintStream();
    Diagnostic.setLogStream(ps.printStream());
    assertEquals(1.0, getLineSolver(true).solveLine(new LinearLine(Double.NEGATIVE_INFINITY, 1.0),  1e-10));
    Diagnostic.setLogStream();
    final String str = ps.toString();
    //System.err.println(str);
    TestUtils.containsAll(str, "starting solveLine", "Infinite case. Defaulting to 1.0", "lo: " + 1.0 + " hi: " + 2.0, "lineValue x: " + 0.0 + " v: -Infinity", "lineValue x: " + 1.0 + " v: -Infinity");
  }

  public void testNaN() {
    assertEquals(0.0, getLineSolver(false).solveLine(new LinearLine(Double.NaN, 1.0), 1e-10));
  }

  private class NastyLine extends Line {
    private final long mConstant;
    NastyLine(final long constant) {
      mConstant = constant;
    }
    @Override
    public double value(final double delta) {
      final long v = Double.doubleToRawLongBits(delta);
      return Double.longBitsToDouble(v ^ (v * mConstant));
    }
  }

  public void testNasty() {
    assertEquals(0.0, getLineSolver(false).solveLine(new NastyLine(3), 1e-10));
    assertEquals(0.0, getLineSolver(false).solveLine(new NastyLine(17), 1e-10));
    assertEquals(0.0, getLineSolver(false).solveLine(new NastyLine(Integer.MAX_VALUE), 1e-10));
  }

  private class ConstantLine extends Line {
    private final double mConstant;
    ConstantLine(final double constant) {
      mConstant = constant;
    }
    @Override
    public double value(final double delta) {
      return mConstant;
    }
  }

  public void testConstant() {
    for (double k = 0; k < 2; ++k) {
      assertEquals(k, getLineSolver(false).solveLine(new ConstantLine(-k), 1e-10));
    }
  }

  public void testInfiniteCL() {
    assertEquals(0.0, getLineSolver(false).solveLine(new ConstantLine(Double.POSITIVE_INFINITY), 1e-10));
    assertEquals(1.0, getLineSolver(false).solveLine(new ConstantLine(Double.NEGATIVE_INFINITY), 1e-10));
  }

  public void testNaNCL() {
    assertEquals(0.0, getLineSolver(false).solveLine(new ConstantLine(Double.NaN), 1e-10));
  }


  private class BizarreLine extends Line {
    private final double mConstant;
    BizarreLine(final double constant) {
      mConstant = constant;
    }
    @Override
    public double value(final double delta) {
      if (delta == 0) {
        return -5;
      }
      if (delta == 1) {
        return 1;
      }
      return mConstant;
    }
  }

  public void testBizarre() {
    assertEquals(0.0, getLineSolver(false).solveLine(new BizarreLine(Double.NaN), 1e-10));
    assertEquals(0.0, getLineSolver(false).solveLine(new BizarreLine(Double.POSITIVE_INFINITY), 1e-10));
    assertEquals(1.0, getLineSolver(false).solveLine(new BizarreLine(Double.NEGATIVE_INFINITY), 1e-10), 1e-10);
    assertEquals(1.0, getLineSolver(false).solveLine(new BizarreLine(-1), 1e-10), 1e-10);
    assertEquals(0.0, getLineSolver(false).solveLine(new BizarreLine(1), 1e-10));
  }

  private class AsymptoteLine extends Line {
    private final double mConstant;
    AsymptoteLine(final double constant) {
      mConstant = constant;
    }
    @Override
    public double value(final double delta) {
      return -1.0 / (delta + mConstant);
    }
  }

  public void testAsymptote() {
    //assertEquals(1.0, LineSolver.solveLine(new AsymptoteLine(0), true, true, System.err), 1e-10);
    assertEquals(1.0, getLineSolver(false).solveLine(new AsymptoteLine(0), 1e-10), 1e-10);
    assertEquals(1.0, getLineSolver(false).solveLine(new AsymptoteLine(1), 1e-10), 1e-10);
  }

  private class NaNPointLine extends Line {
    private final double mPos;
    NaNPointLine(final double pos) {
      mPos = pos;
    }
    @Override
    public double value(final double delta) {
      if (Math.abs(delta - mPos) < 0.00000001) {
        return Double.NaN;
      }
      return delta - 5;
    }
  }

  public void testHitANaN() {
    final MemoryPrintStream ps = new MemoryPrintStream();
    Diagnostic.setLogStream(ps.printStream());
    assertEquals(2.0, getLineSolver(true).solveLine(new NaNPointLine(2), 1e-10), 1e-5);
    Diagnostic.setLogStream();
    final String str = ps.toString();
    assertTrue(str, str.contains("NaN when evaluating hi:"));
  }


  public void testEvalZero() {
    assertEquals(4.0, getLineSolver(false).solveLine(new ZeroLine(), 1e-10), 1e-10);
  }
}
