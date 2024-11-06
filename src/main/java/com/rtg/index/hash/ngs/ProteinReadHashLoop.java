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
package com.rtg.index.hash.ngs;

import java.io.IOException;

import com.rtg.index.hash.ngs.protein.ProteinMask;

/**
 *
 */
public class ProteinReadHashLoop extends ProteinIncrementalHashLoop {

  /**
   * Constructs an instance and as a side effect makes calls to
   * <code>hashCall</code>.
   * @param windowSize number of codes to be included in a window used to make a hash.
   * @param stepSize number of steps to take per hash call
   * @param function used to construct the hash from the codes.
   */
  public ProteinReadHashLoop(int windowSize, int stepSize, ProteinMask function) {
    super(windowSize, stepSize, function, false);
  }

  @Override
  public void hashCall(int internalId, int stepPosition) throws IOException {
    mProteinMask.readAll(internalId, false);
  }

  @Override
  public void hashCallBidirectional(long hashForward, long hashReverse, int stepPosition, int internalId) {
    throw new UnsupportedOperationException("Not supported.");
  }


}
