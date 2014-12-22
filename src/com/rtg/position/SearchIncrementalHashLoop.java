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

import com.rtg.index.Finder;
import com.rtg.index.Index;
import com.rtg.index.hash.ExactHashFunction;
import com.rtg.index.hash.IncrementalHashLoop;
import com.rtg.mode.Frame;
import com.rtg.position.output.PositionOutput;

/**
 */
public class SearchIncrementalHashLoop extends IncrementalHashLoop {

  protected final Finder mHit;

  protected final Finder mHitReverse;

  protected final PositionOutput mOutput;

  protected final PositionOutput mOutputReverse;

  protected final Index mIndex;

  /**
   *
   * @param stepSize step size
   * @param function hash function.
   * @param outputVars variables for output parameters
   * @param outputVarsReverse variables for output parameters (for reverse case)
   * @param index to be updated.
   * @param dualMode whether to process forward and reverse simultaneously or not
   */
  public SearchIncrementalHashLoop(int stepSize, final ExactHashFunction function,
                                   final FinderPositionOutput outputVars, final FinderPositionOutput outputVarsReverse, final Index index, final boolean dualMode)
  {
    super(stepSize, function, dualMode);
    mHit = outputVars.getFinder();
    mOutput = outputVars.getPositionOutput();
    mIndex = index;
    if (outputVarsReverse != null) {
      mOutputReverse = outputVarsReverse.getPositionOutput();
      mHitReverse = outputVarsReverse.getFinder();
    } else {
      mOutputReverse = null;
      mHitReverse = null;
    }
    assert mOutputReverse != mOutput;
  }

  @Override
  public void hashCall(final long hash, final int internalId, final int stepPosition) throws IOException {
    //System.err.println("hashCall hash=" + hash + " id=" + internalId);
    mOutput.setPosition(stepPosition);
    mIndex.search(hash, mHit);
    mOutput.endPosition();
  }

  @Override
  public void hashCallBidirectional(final long hashFoward, final long hashReverse, final int stepPosition, final int internalId) throws IOException {
    hashCall(hashFoward, internalId, stepPosition);
    mOutputReverse.setPosition(stepPosition);
    mIndex.search(hashReverse, mHitReverse);
    mOutputReverse.endPosition();
  }

  @Override
  public void nextSeq(final int seqId, final int length) {
    mOutput.nextSequence(seqId, length);
    if (mOutputReverse != null) {
      mOutputReverse.nextSequence(seqId, length);
    }
  }

  @Override
  public void next(final long seq, final Frame frame) {
    mOutput.nextQuery(frame, (int) seq);
    if (mOutputReverse != null) {
      mOutputReverse.nextQuery(frame.getReverse(), (int) seq);
    }
  }

  @Override
  public void end() throws IOException {
    mOutput.endQuery();
    if (mOutputReverse != null) {
      mOutputReverse.endQuery();
    }
  }

  @Override
  public void endSequence() throws IOException {
    mOutput.endQuerySequence();
    if (mOutputReverse != null) {
      mOutputReverse.endQuerySequence();
    }
  }

  @Override
  public void endAll() throws IOException {
    mOutput.endAll();
    if (mOutputReverse != null) {
      mOutputReverse.endAll();
    }
  }
}
