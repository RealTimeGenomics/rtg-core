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

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.util.array.objectindex.ObjectCreate;
import com.rtg.util.array.objectindex.ObjectIndex;

/**
 * <code>topn</code> implementation for protein
 */
@TestClass("com.rtg.protein.TopNProteinOutputProcessorTest")
public final class TopNProteinImplementation {

  private static final int NUMBER_OF_THREAD_LOCKS = 1 << 16;
  private static final int THREAD_LOCK_MASK = NUMBER_OF_THREAD_LOCKS - 1;

  /* Array of locks for multiple threads */
  private final Object[] mThreadLocks;

  private final short[] mEdgeAlignScore;
  private final byte[] mEdgeScoreCount;
  private final byte[] mResultCount;
  private final ObjectIndex<ProteinAlignmentResult> mResults;
  //private final ProteinAlignmentResult[] mResults;

  private final int mN;

  /**
   * Create a new {@link TopNProteinImplementation}
   * @param topn number of results per read
   * @param numberSequences number of reads
   */
  public TopNProteinImplementation(int topn, int numberSequences) {
    mN = topn;
    final long totalRecords = mN * (long) numberSequences;
    mEdgeAlignScore = new short[numberSequences];
    mResultCount = new byte[numberSequences];
    mEdgeScoreCount = new byte[numberSequences];
    mResults = ObjectCreate.createIndex(totalRecords);
    //mResults = new ProteinAlignmentResult[(int) totalRecords];
    for (int l = 0; l < numberSequences; ++l) {
      mEdgeAlignScore[l] = Short.MAX_VALUE;
    }
    mThreadLocks = new Object[NUMBER_OF_THREAD_LOCKS];
    for (int i = 0; i < mThreadLocks.length; ++i) {
      mThreadLocks[i] = new Object();
    }
  }

  boolean contains(int templateId, int templateStart, int readId, int readAndFrame) {
    for (long i = (long) readId * mN; i < (long) readId * mN + resultCount(readId); ++i) {
      if (mResults.get(i).equals(templateId, templateStart, readAndFrame)) {
        return true;
      }
    }
    return false;
  }

  void insertResult(ProteinAlignmentResult res) {
    final int readId = res.readId();
    synchronized (mThreadLocks[synchId(readId)]) {
      final int currentEdgeScore = mEdgeAlignScore[readId];
      if (res.alignmentScore() < currentEdgeScore) {
        insert(res);
       //insertion
      } else if (res.alignmentScore() == currentEdgeScore) {
        if (mEdgeScoreCount[readId] != -1) {
          mEdgeScoreCount[readId]++; //handle overflow
        }
        if (resultCount(readId) < mN) {
        // insert
          mResults.set((long) readId * mN + resultCount(readId), res);
          mResultCount[readId]++;
        }
      }
    }
  }

  private void insert(ProteinAlignmentResult res) {
    final int readId = res.readId();
    final int alignmentScore = res.alignmentScore();
    final long index = findIndexForInsertion(readId, alignmentScore);
    final long length = mN - (index - (long) readId * mN) - 1;
    updateEdgeCase(readId, alignmentScore);
    for (long i = length - 1; i >= 0; --i) {
      mResults.set(index + i + 1, mResults.get(index + i));
    }
    //System.arraycopy(mResults, index, mResults, index + 1, length);
    mResults.set(index, res);
    mResultCount[readId] = (byte) (resultCount(readId) < mN ? resultCount(readId) + 1 : mN);
  }

  private void updateEdgeCase(int readId, int alignmentScore) {
    final long prevIndex = (long) readId * mN + mN - 2;
    if (!sameScore(mResults.get((long) readId * mN + mN - 1), mResults.get(prevIndex))) {
      final int alignmentScore2 = mResults.get(prevIndex).alignmentScore();
      mEdgeAlignScore[readId] = (short) Math.max(alignmentScore, alignmentScore2);
      if (alignmentScore > alignmentScore2) {
        mEdgeScoreCount[readId] = 1;
      } else {
        final long stop = (long) readId * mN;
        byte count = (byte) (alignmentScore == alignmentScore2 ? 1 : 0);
        for (long i = prevIndex; i >= stop && mResults.get(i).alignmentScore() == alignmentScore2; --i) {
            ++count;
        }
        mEdgeScoreCount[readId] = count;
      }
    }
  }

  private boolean sameScore(ProteinAlignmentResult res1, ProteinAlignmentResult res2) {
    if (res1 == res2) {
      return true;
    }
    if (res1 == null || res2 == null) {
      return false;
    }
    return res1.alignmentScore() == res2.alignmentScore();
  }

  private long findIndexForInsertion(int readId, int alignmentScore) {
    final long indexZero = (long) readId * mN;
    for (long index = (long) readId * mN + resultCount(readId) - 1; index >= indexZero; --index) {
      if (mResults.get(index).alignmentScore() < alignmentScore) {
        return index + 1;
      }
    }
    return indexZero;
  }

  int numResults() {
    return mResultCount.length;
  }

  int resultCount(int readId) {
    return mResultCount[readId] & 0xFF;
  }

  int edgeScore(int readId) {
    return mEdgeAlignScore[readId];
  }

  int edgeScoreCount(int readId) {
    return mEdgeScoreCount[readId] & 0xFF;
  }

  ProteinAlignmentResult result(long id) {
    return mResults.get(id);
  }
  /**
   * Given a read identifier compute the object to synchronize on.
   */
  private int synchId(final int readId) {
    return readId & THREAD_LOCK_MASK;
  }
}
