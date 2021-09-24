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

import com.rtg.metagenomics.matrix.Vector;

/**
 * L values and their derivative along one direction when solving in <code>Species</code>
 */
class LLine extends Line {
  private final SpeciesLine mDeriv;
  private final Vector mR;
  private final Vector mDelta;
  private final BlockInfo mBlockInfo;


  LLine(final Vector r, final Vector delta, final BlockInfo blockInfo) {
    mDeriv = new SpeciesLine(r, delta, blockInfo);
    mR = r;
    mDelta = delta;
    mBlockInfo = blockInfo;
  }

  @Override
  public double value(final double delta) {
    final Vector rp = Species.incrS(mBlockInfo, mR, mDelta, delta);
    return Species.ll(rp, mBlockInfo);
  }

  @Override
  public double[] values(double delta) {
    final double v0 = value(delta);
    final double v1 = mDeriv.value(delta);
    return new double[] {v0, v1};
  }

  @Override
  public int derivativeOrder() {
    return 1;
  }
}
