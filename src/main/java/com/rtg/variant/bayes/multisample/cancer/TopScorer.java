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
package com.rtg.variant.bayes.multisample.cancer;

/**
 * Remember top scoring items.  Simple array backed implementation.
 */
class TopScorer implements Scorer {

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

  @Override
  public void add(final Double score, final int refCount, final int altCount) {
    if (score != null) {
      final double s = score;
      int pos = mLast;
      while (--pos >= 0 && s > mScores[pos]) {
        // do nothing
      }
      ++pos;
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

  @Override
  public int size() {
    return mLast;
  }

  int getRefCount(final int n) {
    return mRefCount[n];
  }

  int getAltCount(final int n) {
    return mAltCount[n];
  }

  private long sum(final int[] array) {
    long s = 0;
    for (final int v : array) {
      s += v;
    }
    return s;
  }

  @Override
  public long getTotalRefCount() {
    return sum(mRefCount);
  }

  @Override
  public long getTotalAltCount() {
    return sum(mAltCount);
  }
}
