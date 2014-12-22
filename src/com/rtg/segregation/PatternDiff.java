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
 * Records the difference between two patterns, in particular if it can be explained by a change of one child.
 */
class PatternDiff extends IntegralAbstract {
  private final int mChild;
  private final boolean mNoDifference;
  private final boolean mNoExplanation;
  private final boolean mFlipped;

  PatternDiff(boolean noDifference, boolean noExplanation, int child, boolean flipped) {
    mChild = child;
    mNoDifference = noDifference;
    mNoExplanation = noExplanation;
    mFlipped = !mNoExplanation && flipped;
    assert integrity();
  }

  int child() {
    if (noDifference() || noExplantion()) {
      throw new RuntimeException("No child defined");
    }
    return mChild;
  }

  boolean noDifference() {
    return mNoDifference;
  }

  boolean noExplantion() {
    return mNoExplanation;
  }

  boolean isValid() {
    return mChild >= 0;
  }

  boolean flipped() {
    if (noExplantion()) {
      throw new RuntimeException("Flip not defined");
    }
    return mFlipped;
  }

  @Override
  public boolean integrity() {
    if (mNoDifference || mNoExplanation) {
      Exam.assertEquals(-1, mChild);
    }
    if (mNoExplanation) {
      Exam.assertFalse(mFlipped);
    }
    return true;
  }
}
