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
      ++mCount;
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
    final long index = mIndex.first(hash);
    if (index < 0) {
      throw new RuntimeException();
    }
    return index;
  }

  @Override
  public boolean contains(Kmer k) {
    final long[] hash = KmerHashA.kmerToHashMin(k);
    if (!mIndex.contains(hash)) {
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
