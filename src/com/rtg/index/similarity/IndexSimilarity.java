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

import java.io.IOException;

import com.rtg.index.Finder;
import com.rtg.index.IndexCompressed;
import com.rtg.index.params.CreateParams;

/**
 */
public final class IndexSimilarity extends IndexCompressed {

  private SimilaritySorter mSimilaritySorter;
  private final boolean mSingleton;

  /**
   * @param indexParams parameters used to size and initialize index.
   * @param threshold maximum number of times a hash can occur before it is ignored.
   * @param proportionalThreshold Whether the frequency threshold should be calculated from index data rather than as a parameter.
   * @param maxThreshold when using proportional threshold don't exceed this repeat frequency
   * @param minThreshold when using proportional threshold don't go below this repeat frequency
   * @param singleton if true then count each hash only once ( as opposed to the product of the number of times it occurs in each genome).
   * @param threads  number of threads appropriate for parallel execution.
   */
  public IndexSimilarity(final CreateParams indexParams, final Integer threshold, final boolean proportionalThreshold, int maxThreshold, int minThreshold, final boolean singleton, final int threads) {
    super(indexParams, threshold, proportionalThreshold, maxThreshold, minThreshold, threads);
    mSingleton = singleton;
  }

  @Override
  public void search(long hash, Finder finder) throws IOException, IllegalStateException {
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
    for (long p = 0; p < mInitialPositionLength - 2; p++) {
      final long hi = mInitialPosition.get(p + 1);
      for (long i = lo; i < hi;) {
      final long hash = mHash.get(i);
      final int seq = (int) mValue.get(i);
      //System.err.println(" i=" + i + " seq=" + seq + " hash=" + hash);
      mSimilaritySorter.add(seq);
      i++;
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
