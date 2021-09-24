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
import com.rtg.index.hash.HashFunction;
import com.rtg.index.hash.ResetHashLoop;
import com.rtg.launcher.HashingRegion;
import com.rtg.launcher.ISequenceParams;
import com.rtg.util.array.WrappedIntArray;

/**
 */
public class BuildResetHashLoop extends ResetHashLoop {

  private final Add mIndex;

  private final long mMxs;

  private int[] mReadLengths = null;

  private final boolean mPairMode;

  private boolean mSecondInPair;


  /**
   * @param windowSize window size
   * @param stepSize step size
   * @param function hash function.
   * @param index to be updated.
   * @param mxs maximum length of sequence.
   */
  public BuildResetHashLoop(final int windowSize, int stepSize, final HashFunction function, final Add index, final long mxs) {
    super(windowSize, stepSize, function, false);
    mIndex = index;
    mMxs = mxs;
    mPairMode = false;
  }

  /**
   * @param windowSize window size
   * @param stepSize step size
   * @param function hash function.
   * @param index to be updated.
   * @param mxs maximum length of sequence.
   * @param pairMode if paired data
   */
  public BuildResetHashLoop(final int windowSize, int stepSize, final HashFunction function, final Add index, final long mxs, final boolean pairMode) {
    super(windowSize, stepSize, function, false);
    mIndex = index;
    mMxs = mxs;
    mPairMode = pairMode;
  }

  @Override
  public long execLoop(ISequenceParams params, final byte[] byteBuffer) throws IOException {
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
    final long length = end - start;
    assert length >= 0;
    if (length > Integer.MAX_VALUE) {
      throw new IllegalArgumentException("too many reads:" + length);
    }
    if (mReadLengths == null) {
      mReadLengths = new int[mPairMode ? (int) length * 2 : (int) length];
    }
    return super.execLoop(params, byteBuffer);
  }

  @Override
  public void nextSeq(final int seqId, final int length) {
    mReadLengths[mPairMode ? (mSecondInPair ? (seqId << 1) + 1 : seqId << 1) : seqId] = length;
  }

  @Override
  public void hashCall(final long hash, final int internalId, final int stepPosition) {
    assert stepPosition >= 0 && stepPosition < mMxs;
    //pairMode translation only works for unidirectional
    final int finId = mPairMode ? (mSecondInPair ? (internalId << 1) + 1 : internalId << 1) : internalId;
    mIndex.add(hash, (long) finId * mMxs + stepPosition);
  }

  @Override
  public WrappedIntArray readLengths() {
    if (mReadLengths == null) {
      throw new IllegalStateException("Read lengths not yet constructed");
    }
    return new WrappedIntArray(mReadLengths);
  }

  @Override
  public void hashCallBidirectional(final long hashForward, final long hashReverse, final int stepPosition, final int internalId) {
    throw new UnsupportedOperationException("Not supported.");
  }

  @Override
  public void setSide(final boolean second) {
    mSecondInPair = second;
  }

}
