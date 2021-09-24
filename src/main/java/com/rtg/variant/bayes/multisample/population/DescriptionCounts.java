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
package com.rtg.variant.bayes.multisample.population;

/**
 * Contains counts of descriptions, useful for computing priors.
 */
public class DescriptionCounts {

  private final int[] mCounts;
  private final int mReferenceIndex;
  private int mTotalCount;

  /**
   * @param size the number of descriptions
   * @param referenceIndex index of description which represents the reference
   */
  public DescriptionCounts(int size, int referenceIndex) {
    mCounts = new int[size];
    mReferenceIndex = referenceIndex;
  }

  void increment(int index, int count) {
    assert index < mCounts.length;
    mCounts[index] += count;
    mTotalCount += count;
  }

  int getReferenceIndex() {
    return mReferenceIndex;
  }

  int getCount(int index) {
    return mCounts[index];
  }

  int[] getCounts() {
    return mCounts;
  }

  int getTotalCount() {
    return mTotalCount;
  }

  public int getSize() {
    return mCounts.length;
  }
}
