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

/**
 * Values along one direction when solving in <code>Species</code>
 */
class SpeciesLineLinearDeriv extends Line {
  private final double mL;
  private final double[] mC;
  private final Frag[] mFrags;


  SpeciesLineLinearDeriv(final Vector r, final Vector delta, final BlockInfo blockInfo) {
    final Vector li = MatrixUtils.pointProduct(blockInfo.getGenomeLengthsVector(), delta);
    mL = MatrixUtils.trace(li);
    //System.err.println("mL=" + mL);
    mFrags = blockInfo.getFrags();
    mC = new double[mFrags.length];
    for (int i = 0; i < mFrags.length; ++i) {
      final Frag frag = mFrags[i];
      final double d = frag.sum(delta);
      final double m = frag.sum(r);
      final double c = d == 0.0 ? 0.0 : m / d;
      mC[i] = c;
      //System.err.println("i=" + i + " d=" + d + " m=" + m + " c=" + c);
    }
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
    double v0 = 0.0;
    double v1 = 0.0;
    for (int i = 0; i < mFrags.length; ++i) {
      final double c = mC[i];
      final int mult = mFrags[i].multiplicity();
      final double t = 1.0 / (c + delta);
      final double tt = mult * t;
      v0 += tt;
      v1 += tt * t;
    }
    final double res0 = mL - v0;
    //System.err.println("delta=" + delta + " v0=" + v0 + "v1=" + v1 + " res=" + res);
    return new double[] {res0, v1};
  }
}
