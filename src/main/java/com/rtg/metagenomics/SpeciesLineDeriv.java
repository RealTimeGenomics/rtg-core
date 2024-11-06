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

import com.rtg.metagenomics.matrix.MatrixUtils;
import com.rtg.metagenomics.matrix.Vector;
import com.rtg.util.diagnostic.Diagnostic;

/**
 * Values and first derivative along one direction when solving in <code>Species</code>
 */
class SpeciesLineDeriv extends Line {
  private final Vector mR;
  final Vector mDelta;
  private final Vector mRD;
  private final Vector mRD2;
  private final Vector mLRD;
  private final Vector mLRD2;
  private final Frag[] mFrags;
  private final int mN;


  SpeciesLineDeriv(final Vector r, final Vector delta, final BlockInfo blockInfo) {
    mR = r;
    mDelta = delta;
    mRD = MatrixUtils.pointProduct(r, delta);
    mRD2 = MatrixUtils.pointProduct(mRD, delta);
    mLRD = MatrixUtils.pointProduct(blockInfo.getGenomeLengthsVector(), mRD);
    mLRD2 = MatrixUtils.pointProduct(mLRD, delta);
    mFrags = blockInfo.getFrags();
    mN = blockInfo.getN();
  }

  @Override
  public double value(double delta) {
    return values(delta)[0];
  }

  @Override
  public int derivativeOrder() {
    return 1;
  }

  @Override
  public double[] values(final double delta) {
    final Vector e = new Vector(mN);
    for (int i = 0; i < mN; ++i) {
      e.set(i, Math.exp(mDelta.get(i) * delta));
    }
    final double v0 = MatrixUtils.multiply(e, mLRD);
    final double v20 = MatrixUtils.multiply(e, mLRD2);
    if (BlockInfo.VERY_VERBOSE) {
      Diagnostic.developerLog("v0=" + v0);
    }

    final Vector re = MatrixUtils.pointProduct(mR, e);
    final Vector red = MatrixUtils.pointProduct(mRD, e);
    final Vector red2 = MatrixUtils.pointProduct(mRD2, e);
    //System.err.println("v0=" + v0 + " re=" + re + " red=" + red);
    double v1 = 0.0;
    double v22 = 0.0;
    double v21 = 0.0;
    for (final Frag frag : mFrags) {
      final double a = frag.sum(red);
      final double a2 = frag.sum(red2);
      final double b = frag.sum(re);
      final double h = a / b;
      //System.err.println("h=" + h + " a=" + a + " b=" + b);
      final int mult = frag.multiplicity();
      v1 -= h * mult;
      v22 += h * h * mult;
      final double h2 = a2 / b;
      v21 -= h2  * mult;
    }
    if (BlockInfo.VERY_VERBOSE) {
      Diagnostic.developerLog("v1=" + v1);
    }
    final double res0 = v0 + v1;
    final double res1 = v20 + v21 + v22;
    //System.err.println("delta=" + delta + " v0=" + v0 + " v1=" + v1 + " res0=" + res0);
    //System.err.println("v20=" + v20 + " v21=" + v21 + " v22=" + v22 + " res1=" + res1);
    return new double[] {res0, res1};
  }
}
