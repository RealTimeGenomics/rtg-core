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
 * Values along one direction when solving in <code>Species</code>
 */
class SpeciesLine extends Line {
  private final Vector mR;
  final Vector mDelta;
  private final Vector mRD;
  private final Vector mLRD;
  private final Frag[] mFrags;
  private final int mN;


  SpeciesLine(final Vector r, final Vector delta, final BlockInfo blockInfo) {
    mR = r;
    mDelta = delta;
    mRD = MatrixUtils.pointProduct(r, delta);
    mLRD = MatrixUtils.pointProduct(blockInfo.getGenomeLengthsVector(), mRD);
    mFrags = blockInfo.getFrags();
    mN = blockInfo.getN();
  }

  @Override
  public double value(final double delta) {
    final Vector e = new Vector(mN);
    for (int i = 0; i < mN; ++i) {
      e.set(i, Math.exp(mDelta.get(i) * delta));
    }
    final double v0 = MatrixUtils.multiply(e, mLRD);
    if (BlockInfo.VERY_VERBOSE) {
      Diagnostic.developerLog("v0=" + v0);
    }

    final Vector re = MatrixUtils.pointProduct(mR, e);
    final Vector red = MatrixUtils.pointProduct(mRD, e);
    //System.err.println("v0=" + v0 + " re=" + re + " red=" + red);
    double v1 = 0.0;
    for (final Frag frag : mFrags) {
      final double a = frag.sum(red);
      final double b = frag.sum(re);
      final double h = a / b;
      //System.err.println("h=" + h + " a=" + a + " b=" + b);
      v1 += h * frag.multiplicity();
    }
    if (BlockInfo.VERY_VERBOSE) {
      Diagnostic.developerLog("v1=" + v1);
    }
    return v0 - v1;
  }
}
