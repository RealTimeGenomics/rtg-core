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
    for (int i = 0; i < mFrags.length; i++) {
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
    for (int i = 0; i < mFrags.length; i++) {
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
