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
    for (int i = mLength; i < newLength; i++) {
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
    for (int i = 0; i < mLength; i++) {
      final int covCnt = mCoverageCount[i];
      for (int j = 0; j <= i; j++) {
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
    for (int i = 0; i < mLength; i++) {
      Exam.assertTrue(mCoverageCount[i] >= 0);
      Exam.assertEquals(i + 1, mBinCount[i].length);
      int count = 0;
      for (int j = 0; j < mBinCount[i].length; j++) {
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
