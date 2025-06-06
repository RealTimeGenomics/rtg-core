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
                                   final FinderPositionOutput outputVars, final FinderPositionOutput outputVarsReverse, final Index index, final boolean dualMode) {
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
  public void hashCallBidirectional(final long hashForward, final long hashReverse, final int stepPosition, final int internalId) throws IOException {
    hashCall(hashForward, internalId, stepPosition);
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
