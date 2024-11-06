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
