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

package com.rtg.variant.bayes.multisample.statistics;

import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 */
class AlleleBalanceTable extends IntegralAbstract {

  private static final String SEP = " ";

  private int mLength = 0;

  private int[] mCoverageCount = new int[mLength];

  private int[][] mBinCount = new int[mLength][];


  /**
   * @param refCount count of number of bases equal to the reference.
   * @param coverage total number of bases.
   */
  public synchronized void update(final int refCount, final int coverage) {
    //System.err.println("update: refCount=" + refCount + " coverage=" + coverage);
    assert refCount <= coverage;
    if (coverage >= mLength) {
      rescale(coverage + 1);
    }
    mCoverageCount[coverage]++;
    mBinCount[coverage][refCount]++;
  }

  private void rescale(int length) {
    int newLength = mLength;
    while (newLength < length) {
      newLength = 2 * newLength + 1;
    }
    final int[] newCoverageCount = new int[newLength];
    System.arraycopy(mCoverageCount, 0, newCoverageCount, 0, mLength);

    final int[][] newBinCount = new int[newLength][];
    System.arraycopy(mBinCount, 0, newBinCount, 0, mLength);
    for (int i = mLength; i < newLength; ++i) {
      newBinCount[i] = new int[i + 1];
    }
    mCoverageCount = newCoverageCount;
    mBinCount = newBinCount;
    mLength = newLength;
  }

  /**
   * Useful only for testing.
   * @return internal length.
   */
  synchronized int length() {
    return mLength;
  }

  @Override
  public synchronized void toString(StringBuilder sb) {
    sb.append("#Coverage" + SEP + "R" + SEP + "X" + SEP + "Count" + SEP + "Total").append(LS);
    for (int i = 0; i < mLength; ++i) {
      final int covCnt = mCoverageCount[i];
      for (int j = 0; j <= i; ++j) {
        final int binCnt = mBinCount[i][j];
        if (binCnt == 0) {
          continue;
        }
        sb.append(i).append(SEP).append(j).append(SEP).append(i - j).append(SEP).append(binCnt).append(SEP).append(covCnt).append(LS);
      }
    }
  }


  @Override
  public synchronized boolean globalIntegrity() {
    integrity();
    for (int i = 0; i < mLength; ++i) {
      Exam.assertTrue(mCoverageCount[i] >= 0);
      Exam.assertEquals(i + 1, mBinCount[i].length);
      int count = 0;
      for (int j = 0; j < mBinCount[i].length; ++j) {
        count += mBinCount[i][j];
      }
      Exam.assertEquals(mCoverageCount[i], count);
    }
    return true;
  }


  @Override
  public synchronized boolean integrity() {
    Exam.assertEquals(mLength, mCoverageCount.length);
    Exam.assertEquals(mLength, mBinCount.length);
    return true;
  }

}
