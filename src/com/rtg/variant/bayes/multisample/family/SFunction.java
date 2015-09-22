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
    for (int i = 0; i < mLength; i++) {
      mForwardS[i + 1] = mForwardS[i].product(mChildS[i]); //TODO use sum if recessive
    }

    mReverseS = new WeightedLattice[mLength + 1];
    mReverseS[mLength] = DefaultWeightedLattice.identity(FastDiseasedFamilyPosterior.SINGLETON, FastDiseasedFamilyPosterior.BIT_SET);
    for (int i = mLength - 1; i >= 0; i--) {
      mReverseS[i] = mReverseS[i + 1].product(mChildS[i]);
    }

    mExcludeS = new WeightedLattice[mLength];
    for (int i = 0; i < mLength; i++) {
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
