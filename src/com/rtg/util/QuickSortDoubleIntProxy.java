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
package com.rtg.util;

/**
 * Quick sort proxy for pair of arrays
 */
public class QuickSortDoubleIntProxy implements QuickSort.SortProxy {
  private final double[] mVals;
  private final int[] mPairs;

  /**
   * Constructor
   * @param valArray values to sort on
   * @param pairArray corresponding array.
   */
  public QuickSortDoubleIntProxy(double[] valArray, int[] pairArray) {
    mVals = valArray;
    mPairs = pairArray;
  }

  @Override
  public int compare(long index1, long index2) {
    return Double.compare(mVals[(int) index2], mVals[(int) index1]);
  }

  @Override
  public long length() {
    return mVals.length;
  }

  @Override
  public void swap(long index1, long index2) {
    final double t = mVals[(int) index1];
    mVals[(int) index1] = mVals[(int) index2];
    mVals[(int) index2] = t;
    final int t2 = mPairs[(int) index1];
    mPairs[(int) index1] = mPairs[(int) index2];
    mPairs[(int) index2] = t2;
  }

  public double[] getVals() {
    return mVals;
  }

  public int[] getPairArray() {
    return mPairs;
  }
}
