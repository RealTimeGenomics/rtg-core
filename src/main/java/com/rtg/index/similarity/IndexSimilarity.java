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
package com.rtg.index.similarity;

import com.rtg.index.Finder;
import com.rtg.index.IndexCompressed;
import com.rtg.index.IndexFilterMethod;
import com.rtg.index.params.CreateParams;

/**
 */
public final class IndexSimilarity extends IndexCompressed {

  private SimilaritySorter mSimilaritySorter;
  private final boolean mSingleton;

  /**
   * @param indexParams parameters used to size and initialize index.
   * @param filter the hash filter
   * @param singleton if true then count each hash only once ( as opposed to the product of the number of times it occurs in each genome).
   * @param threads  number of threads appropriate for parallel execution.
   */
  public IndexSimilarity(final CreateParams indexParams, IndexFilterMethod filter, final boolean singleton, final int threads) {
    super(indexParams, filter, threads);
    mSingleton = singleton;
  }

  @Override
  public void search(long hash, Finder finder) {
    throw new UnsupportedOperationException();
  }

  @Override
  public int count(long hash) {
    throw new UnsupportedOperationException();
  }

  /**
   * Create a self-similarity matrix for the sequences in the index.
   * @param numSequences total number of sequences
   * @return the similarity matrix.
   * @throws IllegalStateException if index has not been frozen.
   */
  public SimilarityMatrix similarity(long numSequences) {
    if (mState != IndexState.FROZEN) {
      throw new IllegalStateException();
    }
    if (mSimilaritySorter == null) {
      mSimilaritySorter = new SimilaritySorter(maxHashCount(), mSingleton);
    }
    final SimilarityMatrix matrix = new SimilarityMatrix(numSequences);
    long lo = 0;
    for (long p = 0; p < mInitialPositionLength - 2; ++p) {
      final long hi = mInitialPosition.get(p + 1);
      for (long i = lo; i < hi;) {
      final long hash = mHash.get(i);
      final int seq = (int) mValue.get(i);
      //System.err.println(" i=" + i + " seq=" + seq + " hash=" + hash);
      mSimilaritySorter.add(seq);
      ++i;
      if (i >= mNumValues) {
        mSimilaritySorter.similarity(matrix);
        //System.err.println(mSimilaritySorter.toString());
        mSimilaritySorter.reset();
        break;
      }
      if (hash != mHash.get(i)) {
        mSimilaritySorter.similarity(matrix);
        //System.err.println(mSimilaritySorter.toString());
        mSimilaritySorter.reset();
      }
    }
      mSimilaritySorter.similarity(matrix);
      //System.err.println(mSimilaritySorter.toString());
      mSimilaritySorter.reset();
      lo = hi;
    }
    return matrix;
  }
}
