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
import com.rtg.index.hash.HashFunction;
import com.rtg.index.hash.ResetHashLoop;
import com.rtg.mode.Frame;
import com.rtg.position.output.PositionOutput;


/**
 */
public class SearchResetHashLoop extends ResetHashLoop {

  private final PositionOutput mOutput;
  private final PositionOutput mOutputReverse;

  private final Finder mHit;
  private final Finder mHitReverse;

  private final Index mIndex;

  /**
   * @param windowSize window size
   * @param stepSize step size
   * @param function hash function.
   * @param outputVars variables for output parameters
   * @param index to be updated.
   * @param dualMode whether to process forward and reverse simultaneously or not
   */
  public SearchResetHashLoop(final int windowSize, int stepSize, final HashFunction function, final FinderPositionOutput outputVars, final Index index, final boolean dualMode) {
    this(windowSize, stepSize, function, outputVars, null, index, dualMode);
  }

  /**
   * @param windowSize window size
   * @param stepSize step size
   * @param function hash function.
   * @param outputVars variables for output parameters
   * @param outputVarsReverse variables for output parameters (for reverse case)
   * @param index to be updated.
   * @param dualMode whether to process forward and reverse simultaneously or not
   */
  public SearchResetHashLoop(int windowSize, int stepSize, HashFunction function, FinderPositionOutput outputVars, FinderPositionOutput outputVarsReverse, Index index, boolean dualMode) {
    super(windowSize, stepSize, function, dualMode);
    mHit = outputVars.getFinder();
    mOutput = outputVars.getPositionOutput();
    mIndex = index;
    if (outputVarsReverse != null) {
      mHitReverse = outputVarsReverse.getFinder();
      mOutputReverse = outputVarsReverse.getPositionOutput();
    } else {
      mHitReverse = null;
      mOutputReverse = null;
    }
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
    mOutput.setPosition(stepPosition);
    mIndex.search(hashFoward, mHit);
    mOutput.endPosition();
    mOutputReverse.setPosition(stepPosition);
    mIndex.search(hashReverse, mHitReverse);
    mOutputReverse.endPosition();
  }

  @Override
  public void nextSeq(final int seqId, final int length) throws IOException {
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
