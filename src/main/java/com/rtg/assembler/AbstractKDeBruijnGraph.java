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

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.index.IndexCompressed;
import com.rtg.index.IndexExtended;
import com.rtg.index.UnfilteredFilterMethod;
import com.rtg.index.params.CreateParams;
import com.rtg.util.MathUtils;
import com.rtg.util.array.bitindex.BitIndex;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.ErrorType;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.diagnostic.OneShotTimer;
import com.rtg.util.iterators.IteratorHelper;

/**
 */
@TestClass("com.rtg.assembler.LowKDeBruijnGraphTest")
public abstract class AbstractKDeBruijnGraph implements DeBruijnGraph {
  protected final IndexExtended mIndex;
  final int mKmerSize;
  protected int mThreshold = 0;
  private final BitIndex mBitIndex;

  @Override
  public long bytes() {
    return mIndex.bytes() + mBitIndex.bytes();
  }
  AbstractKDeBruijnGraph(KmerIterableFactoryInterface factory, final long size, final int kmerSize) {
    mKmerSize = kmerSize;
    final OneShotTimer init = new OneShotTimer("DeBruijn_initial");
    final IndexExtended initialIndex = buildInitialIndex(factory, size, kmerSize);
    final long initialSize = initialIndex.bytes();
    init.stopLog();
    final OneShotTimer count = new OneShotTimer("DeBruijn_count");
    mIndex = buildCountIndex(initialIndex, kmerSize);
    mBitIndex = new BitIndex(mIndex.numberEntries(), 1);
    count.stopLog();
    Diagnostic.developerLog("bytes initial size=" + initialSize + " final size=" + bytes());
  }

  IndexExtended buildInitialIndex(KmerIterableFactoryInterface factory, final long size, final int kmerSize) {
    final int bits = 2 * kmerSize;
    final CreateParams params = new CreateParams(size, bits, bits, 0, true, false, true, true);
    final IndexExtended initialIndex = new IndexCompressed(params, new UnfilteredFilterMethod(), 1);
    Diagnostic.developerLog(initialIndex.infoString());
    try (KmerIterable iterable = factory.makeIterable()) {
      addAll(initialIndex, iterable);
    } catch (IOException e) {
      throw new NoTalkbackSlimException(e, ErrorType.IO_ERROR, e.getMessage());
    }
    try (KmerIterable iterable = factory.makeIterable()) {
      addAll(initialIndex, iterable);
    } catch (IOException e) {
      throw new NoTalkbackSlimException(e, ErrorType.IO_ERROR, e.getMessage());
    }
    return initialIndex;
  }

  private IndexExtended buildCountIndex(final IndexExtended initialIndex, final int kmerSize) {
    final long numberHashes = initialIndex.numberHashes();
    if (numberHashes == 0) {
      // Avoid trying to compact an empty index
      return initialIndex;
    }
    final int maxHashCount = initialIndex.maxHashCount();
    final int valueBits = MathUtils.ceilPowerOf2Bits(maxHashCount);
    final int bits = 2 * kmerSize;
    final CreateParams params = new CreateParams(numberHashes, bits, bits, valueBits, true, true, true, false);
    final IndexExtended countIndex = new IndexCompressed(params, new UnfilteredFilterMethod(), 1);
    Diagnostic.developerLog(countIndex.infoString());

    transferCounts(initialIndex, countIndex);
    transferCounts(initialIndex, countIndex);
    return countIndex;
  }

  protected abstract void transferCounts(final IndexExtended initialIndex, final IndexExtended countIndex);

  private void addAll(final IndexExtended initialIndex, Iterable<Kmer> it) {
    for (final Kmer k : it) {
      add(initialIndex, k);
    }
    initialIndex.freeze();
  }

  protected abstract void add(final IndexExtended initialIndex, Kmer kmer);

  @Override
  public final int frequency(Kmer k) {
    final long search = find(k);
    return (int) mIndex.getValue(search);
  }

  protected abstract long find(Kmer k);

  @Override
  public void setThreshold(int goodThreshold) {
    //TODO optimize by removing entries &le; threshold
    mThreshold  = goodThreshold;
    //System.err.println("threshold=" + mThreshold);
  }

  @Override
  public void setBuilt(Kmer k, boolean built) {
    final long search = find(k);
    mBitIndex.set(search, built ? 1 : 0);
  }

  @Override
  public final boolean isBuilt(Kmer k) {
    final long search = find(k);
    final long built = mBitIndex.get(search);
    assert built == 0 || built == 1;
    return built == 1;
  }

  protected abstract class LocalIterator extends IteratorHelper<Kmer> {
    final long mMaxIndex = mIndex.numberEntries();
    long mNext = 0;

    @Override
    protected void step() {
      ++mNext;
    }

    @Override
    protected boolean atEnd() {
      return mNext >= mMaxIndex;
    }

    @Override
    protected boolean isOK() {
      return mIndex.getValue(mNext) > mThreshold;
    }
  }

}
