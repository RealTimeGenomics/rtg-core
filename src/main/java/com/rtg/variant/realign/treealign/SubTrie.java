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
    ++mTotalCount;
    final int child = start >= bytes.length ? 0 : bytes[start];
    if (start == end || child < 1 || child > 4) {
      ++mStopCount;
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
    for (int i = 0; i < mChildren.length; ++i) {
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
    for (int i = 0; i < 4; ++i) {
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
      Exam.assertTrue(0.0 < mStopProbability && mStopProbability <= 1.0 && !Double.isNaN(mStopProbability));
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
