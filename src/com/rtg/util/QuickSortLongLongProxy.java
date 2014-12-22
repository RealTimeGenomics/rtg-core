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
 * Quick sort proxy for paired long arrays
 */
public class QuickSortLongLongProxy implements QuickSort.SortProxy {
  protected final long[] mVals;
  protected final long[] mPairs;
  protected final boolean mAscending;

  /**
   * Constructor that sorts values in ascending order
   * @param valArray values to sort on
   * @param pairArray corresponding array.
   */
  public QuickSortLongLongProxy(long[] valArray, long[] pairArray) {
    this(valArray, pairArray, true);
  }

  /**
   * Constructor
   * @param valArray values to sort on
   * @param pairArray corresponding array.
   * @param ascending true if sorting should be in ascending order.
   */
  public QuickSortLongLongProxy(long[] valArray, long[] pairArray, boolean ascending) {
    mVals = valArray;
    mPairs = pairArray;
    mAscending = ascending;
  }

  @Override
  public int compare(long index1, long index2) {
    return mAscending
        ? Long.compare(mVals[(int) index1], mVals[(int) index2])
        : Long.compare(mVals[(int) index2], mVals[(int) index1]);
  }

  @Override
  public long length() {
    return mVals.length;
  }

  @Override
  public void swap(long index1, long index2) {
    final long t = mVals[(int) index1];
    mVals[(int) index1] = mVals[(int) index2];
    mVals[(int) index2] = t;
    final long t2 = mPairs[(int) index1];
    mPairs[(int) index1] = mPairs[(int) index2];
    mPairs[(int) index2] = t2;
  }

}
