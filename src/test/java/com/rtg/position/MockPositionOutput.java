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

import com.rtg.mode.Frame;
import com.rtg.position.output.PositionOutput;
import com.rtg.util.integrity.IntegralAbstract;

/**
 */
public class MockPositionOutput extends IntegralAbstract implements PositionOutput {

  @Override
  public void endAll() {
  }

  @Override
  public void endPosition() {
  }

  @Override
  public void endQuery() {
  }

  @Override
  public void endQuerySequence() {
  }

  @Override
  public void hit(final int seqId, final int posn) {
  }

  @Override
  public boolean integrity() {
    return false;
  }

  @Override
  public void nextQuery(final Frame frame, final int seqId) {
  }

  @Override
  public void nextSequence(int seqId, final int length) {
    // do nothing
  }

  @Override
  public void nextSequence(int seqId, final int length, final int usedLength, final byte[] sequence) {
    // do nothing
  }

  @Override
  public void setPosition(final int position) {
  }

  @Override
  public double score() {
    return 0;
  }

  @Override
  public void toString(final StringBuilder sb) {
  }

  /**
   * Return a clone of this suitable for use in the reverse frame
   * @return the clone
   */
  @Override
  public PositionOutput reverseClone() {
    return this;
  }
}
