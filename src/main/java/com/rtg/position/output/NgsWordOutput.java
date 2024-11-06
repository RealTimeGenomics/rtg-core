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
package com.rtg.position.output;

import java.io.IOException;

import com.rtg.index.hash.ngs.OutputProcessor;

/**
 * Simple position output implementation designed to simply feed word hits to an OutputProcessor.
 * Assumes build on read and reverse complement positions have been appropriately adjusted
 */
public final class NgsWordOutput extends AbstractPositionOutput {

  private final OutputProcessor mProcessor;

  /**
   * @param params parameters
   * @param processor processor to output hits to
   */
  public NgsWordOutput(PositionParams params, OutputProcessor processor) {
    super(params, null);
    mProcessor = processor;
  }

  @Override
  public void hit(int seqId, int posn) throws IOException {
    super.hit(seqId, posn);
    final int tPos = mSearchPosition - posn;
    final String frame = mSearchFrame.isForward() ? "F" : "R";
    mProcessor.process(mSearchSeqId, frame, seqId, tPos, 0, 0);
  }

  @Override
  public void memToString(StringBuilder sb) {
  }

  @Override
  public long bytes() {
    return 0;
  }

  @Override
  public void toString(StringBuilder sb) {
  }

  @Override
  public double score() {
    return 0;
  }

  @Override
  public PositionOutput reverseClone() {
    return new NgsWordOutput(mParams, mProcessor);
  }
}
