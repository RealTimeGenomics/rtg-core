/*
 * Copyright (c) 2015. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
 */
package com.rtg.variant.bayes.multisample.cancer;

/**
 * Remember top scoring items.  Simple array backed implementation.
 * @author Sean A. Irvine
 */
class TopScorer {

  // Expect to insert approx 3e6 records into length 1000.  Hence probability
  // of record in result is 1/3000. On average inserted record will be in
  // middle of array, requiring movement of 500 records.  Therefore, expected
  // number of moves if approx. 500/3000 = 5/30 (i.e. less than 1 per record).

  private final double[] mScores;
  private final int[] mRefCount;
  private final int[] mAltCount;
  private int mLast = 0;

  /**
   * @param length maximum values to stores
   */
  TopScorer(final int length) {
    mScores = new double[length];
    mRefCount = new int[length];
    mAltCount = new int[length];
  }

  void add(final Double score, final int refCount, final int altCount) {
    if (score != null) {
      final double s = score;
      int pos = mLast;
      while (--pos >= 0 && s > mScores[pos]) {
        // do nothing
      }
      pos++;
      if (pos < mScores.length) {
        mLast = Math.min(mLast + 1, mScores.length);
        System.arraycopy(mScores, pos, mScores, pos + 1, mLast - pos - 1);
        System.arraycopy(mRefCount, pos, mRefCount, pos + 1, mLast - pos - 1);
        System.arraycopy(mAltCount, pos, mAltCount, pos + 1, mLast - pos - 1);
        mScores[pos] = s;
        mRefCount[pos] = refCount;
        mAltCount[pos] = altCount;
      }
    }
  }

  int size() {
    return mLast;
  }

  int getRefCount(final int n) {
    return mRefCount[n];
  }

  int getAltCount(final int n) {
    return mAltCount[n];
  }
}
