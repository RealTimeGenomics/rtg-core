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
package com.rtg.ngs;

import java.io.IOException;

import com.rtg.index.hash.ngs.OutputProcessor;
import com.rtg.launcher.HashingRegion;

/**
 * A class that can be constructed via <code>threadClone</code> to provide results
 * clipping. The clipping is based on raw hit positions (rather than
 * post-alignment), so should not be used for OutputProcessors that
 * perform alignment themselves.
 *
 */
public class ClippedOutputProcessor implements OutputProcessor {

  private final OutputProcessor mInner;

  private final HashingRegion mRegion;


  /**
   * Creates a new <code>ClippedOutputProcessor</code> instance.
   *
   * @param inner an <code>OutputProcessor</code> to which non-filtered hits are delegated
   * @param region a <code>ClipRegion</code> value
   */
  public ClippedOutputProcessor(OutputProcessor inner, HashingRegion region) {
    mInner = inner;
    mRegion = region;
  }

  // Implementation of com.rtg.index.hash.ngs.OutputProcessor

  @Override
  public void process(final long templateId, final String frame, final int readId, final int tStart, final int score, final int scoreIndel) throws IOException {
    if (mRegion.isInRange(templateId, tStart)) {
      mInner.process(templateId, frame, readId, tStart, score, scoreIndel);
    }
  }

  @Override
  public void finish() {
    // Noop
  }

  @Override
  public OutputProcessor threadClone(final HashingRegion clipRegion) {
    throw new UnsupportedOperationException();
  }
  @Override
  public void threadFinish() throws IOException {
    mInner.threadFinish();
  }

  @Override
  public void close() {
    // Noop
  }

}
