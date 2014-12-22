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
import com.rtg.metagenomics.matrix.VectorSimple;
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
    final Vector e = new VectorSimple(mN);
    for (int i = 0; i < mN; i++) {
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
