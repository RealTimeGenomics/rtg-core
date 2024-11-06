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

package com.rtg.reader;

import java.io.Closeable;
import java.io.IOException;
import java.util.function.Predicate;

import com.rtg.util.Histogram;

import htsjdk.samtools.util.RuntimeIOException;

/**
 * Chunk writer for {@code FastqPair}
 */
class AsyncFastqPairWriter extends AbstractAsyncChunkWriter<FastqPair> {

  private final FastqWriter mLeft;
  private final FastqWriter mRight;
  private final Predicate<FastqPair> mAccept;
  private final Histogram mLeftLengths;
  private final Histogram mRightLengths;


  AsyncFastqPairWriter(FastqWriter r1, FastqWriter r2) {
    this(r1, r2, o -> true);
  }

  AsyncFastqPairWriter(FastqWriter r1, FastqWriter r2, Predicate<FastqPair> accept) {
    this(r1, r2, accept, null, null);
  }

  AsyncFastqPairWriter(FastqWriter r1, FastqWriter r2, Predicate<FastqPair> accept, Histogram leftHist, Histogram rightHist) {
    super(10000);
    mLeft = r1;
    mRight = r2;
    mAccept = accept;
    mLeftLengths = leftHist;
    mRightLengths = rightHist;
  }

  @Override
  public void write(FastqPair pair) {
    if (mAccept.test(pair)) {
      super.write(pair);
    }
  }

  @Override
  protected String getThreadNamePrefix() {
    return "FastqWriter";
  }

  @Override
  protected void synchronouslyWrite(FastqPair pair) {
    try {
      if (mLeftLengths != null) {
        mLeftLengths.increment(pair.r1().length());
      }
      if (mRightLengths != null) {
        mRightLengths.increment(pair.r2().length());
      }
      synchronouslyWrite(mLeft, pair.r1());
      synchronouslyWrite(mRight, pair.r2());
    } catch (IOException e) {
      throw new RuntimeIOException(e);
    }
  }

  private void synchronouslyWrite(FastqWriter w, FastqSequence read) throws IOException {
    w.write(read);
  }

  @Override
  @SuppressWarnings("try")
  protected void synchronouslyClose() {
    try {
      try (Closeable ignored = mLeft; Closeable ignored2 = mRight) {
        // Just for nice closing behaviour
      }
    } catch (IOException e) {
      throw new RuntimeIOException(e);
    }
  }
}
