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
package com.rtg.index.similarity;

import java.util.Arrays;

import com.rtg.util.StringUtils;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 * Takes a set or sequence ids and sorts them removes duplicates and then
 * updates the similarity matrix. This is a potentially very efficient
 * way of doing this especially when there are a large number of sequence ids
 * with many duplicates.
 */
public final class SimilaritySorter extends IntegralAbstract {

  private final int mLength;

  private final boolean mSingleton;

  private int mCurr;

  private final int[] mSortedIds;

  private boolean mSorted;

  private final long[] mCounts;

  private int mNumNoDupl;

  private final int[] mNoDupl;

  @Override
  public boolean globalIntegrity() {
    integrity();
    final int[] sId = new int[mCurr];
    for (int i = 0; i < mCurr; i++) {
      Exam.assertTrue(mSortedIds[i] >= 0);
      sId[i] = mSortedIds[i];
    }
    if (mSorted) {
      //original array in order
      for (int i = 1; i < mCurr; i++) {
        Exam.assertTrue(mSortedIds[i - 1] <= mSortedIds[i]);
      }
      final int[] sNo = new int[mNumNoDupl];
      for (int i = 0; i < mNumNoDupl; i++) {
        Exam.assertTrue(mNoDupl[i] >= 0);
        sNo[i] = mNoDupl[i];
        //make sure everything in no duplicates is in original
        Exam.assertTrue(Arrays.binarySearch(sId, mNoDupl[i]) >= 0);
      }
      //no duplicates in order
      for (int i = 1; i < mNumNoDupl; i++) {
        Exam.assertTrue(mNoDupl[i - 1] < mNoDupl[i]);
      }
      //make sure every original occurs in the no duplicates
      for (int i = 0; i < mCurr; i++) {
        Exam.assertTrue(String.valueOf(mSortedIds[i]), Arrays.binarySearch(sNo, mSortedIds[i]) >= 0);
      }
    }
    return true;
  }

  @Override
  public boolean integrity() {
    //Assert.assertTrue(mLength > 0);
    Exam.assertTrue(mCurr >= 0 && mCurr <= mLength);
    if (mSorted) {
      Exam.assertTrue(mNumNoDupl <= mCurr);
    } else {
      Exam.assertEquals(-1, mNumNoDupl);
    }
    return true;
  }


  /**
   * @param length maximum number of sequence ids that can be recorded.
   * @param singleton if true then count each hash only once ( as opposed to the product of the number of times it occurs in each genome).
   */
  public SimilaritySorter(final int length, final boolean singleton) {
    mLength = length;
    mSingleton = singleton;
    mSortedIds = new int[mLength];
    mCounts = new long[mLength];
    mNoDupl = new int[mLength];
    reset();
    integrity();
  }

  /**
   * Add a new sequence id.
   * @param seqId sequence id.
   */
  public void add(final int seqId) {
    assert !mSorted;
    mSortedIds[mCurr] = seqId;
    mCurr++;
  }

  /**
   * Reset.
   */
  public void reset() {
    mSorted = false;
    mCurr = 0;
    mNumNoDupl = -1;
  }

  /**
   * Increment the matrix using accumulated values.
   * @param matrix to be updated.
   */
  public void similarity(final SimilarityMatrix matrix) {
    Arrays.sort(mSortedIds, 0, mCurr);
//    System.err.print("[");
//    for (int i = 0; i < mCurr; i++) {
//      System.err.print(mSortedIds[i] + " ");
//    }
//    System.err.println("]");
    //System.err.println(Arrays.toString(mSortedIds));
    //remove duplicates and count them
    int to = -1;
    int last = -1;
    for (int from = 0; from < mCurr; from++) {
      final int fromV = mSortedIds[from];
      if (fromV != last) {
        to++;
        mNoDupl[to] = fromV;
        mCounts[to] = 0;
      }
      mCounts[to]++;
      last = fromV;
    }
    mNumNoDupl = to + 1;

    //put counts into similarity matrix
    for (int j = 0; j < mNumNoDupl; j++) {
      for (int k = j; k < mNumNoDupl; k++) {
        if (mSingleton) {
          matrix.increment(mNoDupl[j], mNoDupl[k]);
        } else {
          matrix.increment(mNoDupl[j], mNoDupl[k], (double) mCounts[j] * mCounts[k]); //TODO add different ways of computing the value to increment by
        }
      }
    }
    mSorted = true;
  }

  @Override
  public void toString(final StringBuilder sb) {
    sb.append("SimilaritySorter ").append(mCurr).append(":").append(mLength).append(" ").append(mSorted ? "sorted" : "unsorted").append(StringUtils.LS);
    for (int i = 0; i < mCurr; i++) {
      sb.append(mSortedIds[i]).append(" ");
    }
    sb.append(StringUtils.LS);
    if (mSorted) {
      sb.append("duplicates removed ").append(mNumNoDupl).append(StringUtils.LS);
      for (int i = 0; i < mNumNoDupl; i++) {
        sb.append(mNoDupl[i]).append(" ");
      }
      sb.append(StringUtils.LS);
      for (int i = 0; i < mNumNoDupl; i++) {
        sb.append(mCounts[i]).append(" ");
      }
      sb.append(StringUtils.LS);
    }
  }
}
