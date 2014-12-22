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

package com.rtg.variant.realign.treealign;

import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 *
 */
class SubTrie extends IntegralAbstract {
  int mTotalCount;
  int mStopCount;
  final SubTrie[] mChildren = new SubTrie[4];

  double mStopProbability = Double.NaN;
  boolean mFrozen = false;

  SubTrie(final int total, final int stop) {
    mTotalCount = total;
    mStopCount = stop;
  }

  SubTrie() {
    // dont need to do anything
  }

  void increment(final byte[] bytes, final int start, final int end) {
    assert !mFrozen;
    mTotalCount++;
    final int child = start >= bytes.length ? 0 : bytes[start];
    if (start == end || child < 1 || child > 4) {
      mStopCount++;
      return;
    }
    Trie.increment(mChildren, child - 1, bytes, start + 1, end);
  }

  @Override
  public boolean equals(Object obj) {
    if (obj == null) {
      return false;
    }
    if (this == obj) {
      return true;
    }
    final SubTrie that = (SubTrie) obj;
    if (this.mTotalCount != that.mTotalCount || this.mStopCount != that.mStopCount) {
      return false;
    }
    for (int i = 0; i < mChildren.length; i++) {
      final SubTrie thisCh = this.mChildren[i];
      final SubTrie thatCh = that.mChildren[i];
      if (thisCh == null && thatCh == null) {
        continue;
      }
      if (!(thisCh != null && thatCh != null)) {
        return false;
      }
      if (!thisCh.equals(thatCh)) {
        return false;
      }
    }
    return true;
  }

  @Override
  public int hashCode() {
    return super.hashCode();
  }

  @Override
  public boolean globalIntegrity() {
    integrity();
    int cnt = mStopCount;
    for (int i = 0; i < 4; i++) {
      final SubTrie child = mChildren[i];
      if (child != null) {
        cnt += child.mTotalCount;
      }
    }
    Exam.assertEquals(mTotalCount, cnt);
    return true;
  }

  @Override
  public boolean integrity() {
    Exam.assertTrue(mStopCount >= 0);
    Exam.assertTrue(mTotalCount >= mStopCount);
    Exam.assertEquals(4, mChildren.length);
    if (!mFrozen) {
      Exam.assertTrue(Double.isNaN(mStopProbability));
    } else {
      Exam.assertTrue("" + mStopProbability, 0.0 < mStopProbability && mStopProbability <= 1.0 && !Double.isNaN(mStopProbability));
    }
    return true;
  }

  /**
   * @param priorStop prior probability of stopping.
   * @param soFar the probability so far down this branch.
   */
  void freeze(final double priorStop, final double soFar) {
    assert !mFrozen;
    final double nullPrior = (1.0 - priorStop) / 4.0;
    double nonNull = mTotalCount - mStopCount;
    double nullPriorTotal = 0.0;
    for (SubTrie aMChildren : mChildren) {
      if (aMChildren == null) {
        nullPriorTotal += nullPrior;
      } else {
        nonNull += nullPrior;
      }
    }
    final double stop = mStopCount + priorStop + nullPriorTotal;
    final double total = mTotalCount + 1.0;
    mStopProbability = soFar * (stop / total);
    assert 0.0 < mStopProbability && mStopProbability <= 1.0 && !Double.isNaN(mStopProbability);
    final double cont0 = soFar - mStopProbability;
    for (final SubTrie child : mChildren) {
      if (child != null) {
        final double childWeight = child.mTotalCount + nullPrior;
        final double childSofar = cont0 * childWeight / nonNull;
        child.freeze(priorStop, childSofar);
      }
    }
    mFrozen = true;
  }

}
