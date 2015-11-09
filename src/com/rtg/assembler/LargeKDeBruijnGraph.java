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

import java.util.Arrays;
import java.util.Iterator;

import com.rtg.index.FinderHashValueExtended;
import com.rtg.index.IndexExtended;

/**
 * Deals with case where 32 &lt; k.
 */
public class LargeKDeBruijnGraph extends AbstractKDeBruijnGraph {

  LargeKDeBruijnGraph(KmerIterableFactoryInterface factory, final long size, final int kmerSize) {
    super(factory, size, kmerSize);
  }

  private static class LocalFinderHashValue implements FinderHashValueExtended {
    private final IndexExtended mCountIndex;

    int mCount = 0;
    boolean mFirst = true;
    long[] mLastHash;

    LocalFinderHashValue(IndexExtended countIndex) {
      super();
      mCountIndex = countIndex;
    }

    @Override
    public void found(long[] hash, long value) {
      if (!Arrays.equals(hash, mLastHash)) {
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
      if (!mFirst) {
        mCountIndex.add(mLastHash, mCount);
      }
    }

  }

  @Override
  protected void transferCounts(final IndexExtended initialIndex, final IndexExtended countIndex) {
    try {
      final LocalFinderHashValue finder = new LocalFinderHashValue(countIndex);
      initialIndex.scanAll(finder);
      finder.atEnd();
      countIndex.freeze();
    } catch (final Exception e) {
      throw new RuntimeException(e);
    }

  }

  @Override
  protected final void add(final IndexExtended initialIndex, Kmer kmer) {
    final long[] hash = KmerHashA.kmerToHashMin(kmer);
    //System.err.println("add kmer=" + kmer);
    //System.err.println("    hash=" + Utils.toBits(hash));
    initialIndex.add(hash, 0L);
  }

  @Override
  protected final long find(Kmer k) {
    final long[] hash = KmerHashA.kmerToHashMin(k);
    final long search = mIndex.contains(hash);
    if (search < 0) {
      throw new RuntimeException();
    }
    return search;
  }

  @Override
  public boolean contains(Kmer k) {
    final long[] hash = KmerHashA.kmerToHashMin(k);
    final long search = mIndex.contains(hash);
    if (search < 0) {
      return false;
    }
    final int frequency = frequency(k);
    return frequency > mThreshold;
  }

  private class LocalIteratorLargek extends LocalIterator {
    @Override
    protected Kmer current() {
      final long[] hash = mIndex.getHashExtended(mNext);
      //System.err.println("get kmer=" + kmer);
      //System.err.println("    hash=" + Utils.toBits(hash));
      return new KmerHashA(hash, mKmerSize);
    }
  }

  @Override
  public Iterator<Kmer> iterator() {
    return new LocalIteratorLargek();
  }

}
