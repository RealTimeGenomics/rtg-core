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
package com.rtg.position;

import java.io.IOException;

import com.rtg.index.Add;
import com.rtg.index.hash.ExactHashFunction;
import com.rtg.index.hash.IncrementalHashLoop;
import com.rtg.launcher.HashingRegion;
import com.rtg.launcher.ISequenceParams;
import com.rtg.util.array.ImmutableIntArray;
import com.rtg.util.array.SingleValueIntArray;
import com.rtg.util.array.WrappedIntArray;

/**
 */
public class BuildIncrementalHashLoop extends IncrementalHashLoop {

  private final Add mIndex;

  private final int mMxs;

  private int[] mReadLengths = null;
  private int mReadLength = -1;
  private int mNumReads = 0;

  private final boolean mPairMode;

  private boolean mSecondInPair;

  /**
   * @param stepSize step size
   * @param function hash function.
   * @param index to be updated.
   * @param mxs maximum length of sequence.
   * @param pairMode if paired data
   */
  public BuildIncrementalHashLoop(int stepSize, final ExactHashFunction function, final Add index, final int mxs, final boolean pairMode) {
    super(stepSize, function, false);
    mIndex = index;
    mMxs = mxs;
    mPairMode = pairMode;
  }

  @Override
  public long execLoop(ISequenceParams params, byte[] byteBuffer)
      throws IOException {
    final HashingRegion region = params.region();
    final long start;
    final long end;
    if (region != HashingRegion.NONE) {
      start = region.getStart();
      end = region.getEnd();
    } else {
      start = 0;
      end = params.reader().numberSequences();
    }
    final long len = end - start;
    assert len >= 0;
    if (len > Integer.MAX_VALUE) {
      throw new IllegalArgumentException("too many reads:" + len);
    }
    if (mReadLengths == null) {
      if (params.reader().maxLength() == params.reader().minLength()) {
        mReadLength = (int) params.reader().maxLength();
        mNumReads = mPairMode ? (int) len * 2 : (int) len;
      } else {
        mReadLengths = new int[mPairMode ? (int) len * 2 : (int) len];
      }
    }
    return super.execLoop(params, byteBuffer);
  }


  @Override
  public void nextSeq(final int seqId, final int length) {
    //System.err.println("nextSeq seqId=" + seqId + " length=" + length);
    if (mReadLengths != null) {
      mReadLengths[mPairMode ? (mSecondInPair ? (seqId << 1) + 1 : seqId << 1) : seqId] = length;
    }
  }

  @Override
  public void hashCall(final long hash, final int internalId, final int stepPosition) {
    //System.err.println("build hashCall hash=" + hash + " id=" + internalId);
    assert stepPosition >= 0 && stepPosition < mMxs;
    //pairMode translation only works for unidirectional
    final int finId = mPairMode ? (mSecondInPair ? (internalId << 1) + 1 : internalId << 1) : internalId;
    final long id = (long) finId * mMxs + stepPosition;
    //System.err.println("build hashCall hash=" + hash + " id=" + internalId + " actual id=" + id + " finId=" + finId + " mxs=" + mMxs + " stepPosition=" + stepPosition);
    mIndex.add(hash, id);
  }

  @Override
  public ImmutableIntArray readLengths() {
    if (mReadLengths == null && mReadLength == -1) {
      throw new IllegalStateException("Read lengths not yet constructed");
    }
    if (mReadLength >= 0) {
      return new SingleValueIntArray(mReadLength, mNumReads);
    } else {
      return new WrappedIntArray(mReadLengths);
    }
  }

  @Override
  public void hashCallBidirectional(final long hashFoward, final long hashReverse, final int stepPosition, final int internalId) {
    throw new UnsupportedOperationException("Not supported.");
  }

  @Override
  public void setSide(final boolean second) {
    mSecondInPair = second;
  }
}
