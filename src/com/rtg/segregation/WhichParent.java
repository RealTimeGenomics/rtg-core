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

package com.rtg.segregation;

import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 * Records if valid single parent and child which explains a cross-over.
 */
class WhichParent extends IntegralAbstract {
  private final PatternDiff mFather;
  private final PatternDiff mMother;
  private final boolean mValid;

  WhichParent(PatternDiff father, PatternDiff mother) {
    mFather = father;
    mMother = mother;
    mValid = father.noDifference() && mother.isValid() || mother.noDifference() && father.isValid();
    assert integrity();
  }

  boolean isValid() {
    return mValid;
  }

  boolean isFather() {
    return mFather.isValid();
  }

  int child() {
    return isFather() ? mFather.child() : mMother.child();
  }

  PatternDiff father() {
    return mFather;
  }

  PatternDiff mother() {
    return mMother;
  }

  @Override
  public boolean integrity() {
    Exam.assertNotNull(mFather);
    Exam.assertNotNull(mMother);
    return true;
  }
}
