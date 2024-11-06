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
