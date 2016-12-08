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
package com.rtg.protein;

import com.rtg.util.array.objectindex.ObjectCreate;
import com.rtg.util.array.objectindex.ObjectIndex;

/**
 * Top equal implementation for proteins
 */
public final class TopEqualProteinImplementation {

  private static final int NUMBER_OF_THREAD_LOCKS = 1 << 16;
  private static final int THREAD_LOCK_MASK = NUMBER_OF_THREAD_LOCKS - 1;

  /* Array of locks for multiple threads */
  private static final Object[] THREAD_LOCKS;

  static {
    THREAD_LOCKS = new Object[NUMBER_OF_THREAD_LOCKS];
    for (int i = 0; i < THREAD_LOCKS.length; ++i) {
      THREAD_LOCKS[i] = new Object();
    }
  }

  private final short[] mBestAlignScore;
  private final byte[] mResultCount;
  private final ObjectIndex<ProteinAlignmentResult> mResults;
  //private final ProteinAlignmentResult[] mResults;
  private final int mTopN;

  /**
   * Construct a {@link TopEqualProteinImplementation}
   * @param topn number of results
   * @param numberSequences number of sequences
   */
  protected TopEqualProteinImplementation(int topn, int numberSequences) {
    final long totalRecords = topn * (long) numberSequences;
    mTopN = topn;
    mBestAlignScore = new short[numberSequences];
    mResultCount = new byte[numberSequences];
    mResults = ObjectCreate.createIndex(totalRecords);
    //mResults = new ProteinAlignmentResult[(int) totalRecords];
    for (int l = 0; l < numberSequences; ++l) {
      mBestAlignScore[l] = Short.MAX_VALUE;
    }
  }

  boolean contains(int templateId, int templateStart, int readId, int readAndFrame) {
    synchronized (THREAD_LOCKS[synchId(readId)]) {
      final int count = resultCount(readId) > mTopN ? mTopN : resultCount(readId);
      for (long i = (long) readId * mTopN; i < (long) readId * mTopN + count; ++i) {
        if (mResults.get(i).equals(templateId, templateStart, readAndFrame)) {
          return true;
        }
      }
      return false;
    }
  }

  void insertResult(ProteinAlignmentResult res) {
    final int readId = res.readId();
    synchronized (THREAD_LOCKS[synchId(readId)]) {
      final short currentBestScore = mBestAlignScore[readId];
      if (res.alignmentScore() < currentBestScore) {
        mBestAlignScore[readId] = (short) res.alignmentScore();
        assert res.alignmentScore() == mBestAlignScore[readId];
        mResultCount[readId] = (byte) 1;
        mResults.set((long) readId * mTopN, res);
      } else if (res.alignmentScore() == currentBestScore) {
        final int currentCount = resultCount(readId);
        if (currentCount != 255) {
          mResultCount[readId]++;
          if (currentCount < mTopN) {
            mResults.set((long) readId * mTopN + currentCount, res);
          }
        }
      }
    }
  }

  int numResults() {
    return mResultCount.length;
  }

  int resultCount(int readId) {
    return mResultCount[readId] & 0xFF;
  }

  int score(final int readId) {
    return mBestAlignScore[readId];
  }

  ProteinAlignmentResult result(long resultId) {
    return mResults.get(resultId);
  }

  /**
   * Given a read identifier compute the object to synchronize on.
   *
   */
  private int synchId(final int readId) {
    return readId & THREAD_LOCK_MASK;
  }
}
