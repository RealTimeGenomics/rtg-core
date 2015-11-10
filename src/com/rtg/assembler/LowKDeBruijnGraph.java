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

package com.rtg.assembler;

import java.io.IOException;
import java.util.Iterator;

import com.rtg.index.FinderHashValue;
import com.rtg.index.Index;
import com.rtg.index.IndexExtended;

/**
 * Deals with case where k &le; 32.
 */
public class LowKDeBruijnGraph extends AbstractKDeBruijnGraph {

  LowKDeBruijnGraph(KmerIterableFactoryInterface factory, final long size, final int kmerSize) {
    super(factory, size, kmerSize);
    assert kmerSize <= 32;
  }

  private static class LocalFinderHashValue implements FinderHashValue {
    private final Index mCountIndex;

    int mCount = 0;
    boolean mFirst = true;
    long mLastHash;

    LocalFinderHashValue(Index countIndex) {
      super();
      mCountIndex = countIndex;
    }

    @Override
    public void found(long hash, long value) {
      if (hash != mLastHash) {
        if (!mFirst) {
          mCountIndex.add(mLastHash, mCount);
        }
        mCount = 0;
        mLastHash = hash;
      }
      mFirst = false;
      mCount++;
    }

    void atEnd() {
      mCountIndex.add(mLastHash, mCount);
    }

  }

  @Override
  protected void transferCounts(final IndexExtended initialIndex, final IndexExtended countIndex) {
    try {
      final LocalFinderHashValue finder = new LocalFinderHashValue(countIndex);
      initialIndex.scan(finder);
      finder.atEnd();
      countIndex.freeze();
    } catch (final IOException e) {
      throw new RuntimeException(e);
    }

  }

  @Override
  protected final void add(final IndexExtended initialIndex, Kmer kmer) {
    final long hash = KmerHash.kmerToHashMin(kmer);
    initialIndex.add(hash, 0L);
  }

  @Override
  protected final long find(Kmer k) {
    final long hash = KmerHash.kmerToHashMin(k);
    final long search = mIndex.first(hash);
    if (search < 0) {
      throw new RuntimeException();
    }
    return search;
  }

  @Override
  public boolean contains(Kmer k) {
    //System.err.println("contains(" + k + ")");
    final long hash = KmerHash.kmerToHashMin(k);
    if (!mIndex.contains(hash)) {
      return false;
    }
    final int frequency = frequency(k);
    //System.err.println("frequency=" + frequency + " threshold=" + mThreshold);
    return frequency > mThreshold;
  }

  private class LocalIteratorLowk extends LocalIterator {
    @Override
    protected Kmer current() {
      final long hash = mIndex.getHash(mNext);
      return new KmerHash(hash, mKmerSize);
    }
  }

  @Override
  public Iterator<Kmer> iterator() {
    return new LocalIteratorLowk();
  }

}
