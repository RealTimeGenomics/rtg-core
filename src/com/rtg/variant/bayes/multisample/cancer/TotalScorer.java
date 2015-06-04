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
 * Remember total of all items.
 */
class TotalScorer implements Scorer {

  private long mRefCount = 0;
  private long mAltCount = 0;
  private int mCount = 0;

  @Override
  public void add(final Double score, final int refCount, final int altCount) {
    if (score != null) {
      mCount++;
      mRefCount += refCount;
      mAltCount += altCount;
    }
  }

  @Override
  public int size() {
    return mCount;
  }

  @Override
  public long getTotalRefCount() {
    return mRefCount;
  }

  @Override
  public long getTotalAltCount() {
    return mAltCount;
  }
}
