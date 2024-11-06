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
 * Wraps around a line solver implementation adding developer logging.
 */
final class LoggingLineSolver implements LineSolver {

  static final double RELATIVE_THRESHOLD = 0.01;

  private final boolean mVeryVerbose;

  private final LineSolver mSolver;

  /**
   * Constructor setting verbose debugging output.
   * @param solver the delegate line solver
   * @param veryVerbose verbose output.
   */
  LoggingLineSolver(LineSolver solver, boolean veryVerbose) {
    if (solver == null) {
      throw new NullPointerException("null solver given.");
    }
    mSolver = solver;
    mVeryVerbose = veryVerbose;
  }

  @Override
  public double solveLine(final Line line, final double relThreshold) {
    final Line vline = new LineVerbose(line, mVeryVerbose);
    final double x = mSolver.solveLine(vline, relThreshold);

    assert x >= 0.0 && AbstractSolver.isFinite(x);
    //assert x == 0.0 || line.value(x) <= 0.0;
    return x;
  }

  private static class LineVerbose extends Line {
    private final Line mProxy;
    private final boolean mVeryVerbose;
    /**
     * @param proxy proxy line
     * @param veryVerbose print very verbose information
     */
    LineVerbose(Line proxy, boolean veryVerbose) {
      super();
      mProxy = proxy;
      mVeryVerbose = veryVerbose;
    }
    @Override
    public double value(double x) {
      final double v = mProxy.value(x);
      if (mVeryVerbose) {
        Diagnostic.developerLog("lineValue x: " + x + " v: " + v);
      }
      return v;
    }
  }
}
