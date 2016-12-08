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
