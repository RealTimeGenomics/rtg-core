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
package com.rtg.variant.bayes.multisample.family;

import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 * Incrementally compute function <code>S</code> as defined in the multi-calling
 * documentation.
 */
//TODO version that uses sum rather than product for recessive rather than dominant disease
public class SFunction extends IntegralAbstract {

  private final WeightedLattice mInit;
  private final WeightedLattice[] mChildS;
  private final WeightedLattice[] mForwardS;
  private final WeightedLattice[] mReverseS;
  private final WeightedLattice[] mExcludeS;
  private final int mLength;

  //TODO use flag or separate version for recessive.
  SFunction(final WeightedLattice init, final WeightedLattice[] children) {
    mInit = init;
    mLength = children.length;
    mChildS = children;

    mForwardS = new WeightedLattice[mLength + 1];
    mForwardS[0] = init;
    for (int i = 0; i < mLength; ++i) {
      mForwardS[i + 1] = mForwardS[i].product(mChildS[i]); //TODO use sum if recessive
    }

    mReverseS = new WeightedLattice[mLength + 1];
    mReverseS[mLength] = DefaultWeightedLattice.identity(FastDiseasedFamilyPosterior.SINGLETON, FastDiseasedFamilyPosterior.BIT_SET);
    for (int i = mLength - 1; i >= 0; --i) {
      mReverseS[i] = mReverseS[i + 1].product(mChildS[i]);
    }

    mExcludeS = new WeightedLattice[mLength];
    for (int i = 0; i < mLength; ++i) {
      mExcludeS[i] = mForwardS[i].product(mReverseS[i + 1]);
    }
    assert globalIntegrity();
  }

  WeightedLattice all() {
    return mForwardS[mLength];
  }

  WeightedLattice excludeChild(final int i) {
    return mExcludeS[i];
  }

  @Override
  public final boolean globalIntegrity() {
    integrity();
    final WeightedLattice f = mForwardS[mLength];
    f.globalIntegrity();
    final WeightedLattice r = mInit.product(mReverseS[0]);
    r.globalIntegrity();
    Exam.assertTrue(f.approxEqual(r, 0.0001));
    return true;
  }

  @Override
  public boolean integrity() {
    Exam.assertEquals(mLength, mChildS.length);
    Exam.assertEquals(mLength + 1, mForwardS.length);
    Exam.assertEquals(mLength + 1, mReverseS.length);
    Exam.assertEquals(mLength, mExcludeS.length);
    return true;
  }

}
