/*
 * Copyright (c) 2016. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
 */

package com.rtg.reader;

import java.io.Closeable;
import java.io.IOException;

import htsjdk.samtools.util.RuntimeIOException;

/**
 * Chunk writer for {@code FastqPair}
 */
class AsyncFastqPairWriter extends AbstractAsyncChunkWriter<FastqPair> {
  private final FastqWriter mLeft;
  private final FastqWriter mRight;

  AsyncFastqPairWriter(FastqWriter r1, FastqWriter r2) {
    super(10000);
    mLeft = r1;
    mRight = r2;
  }

  @Override
  protected String getThreadNamePrefix() {
    return "FastqWriter";
  }

  @Override
  protected void synchronouslyWrite(FastqPair pair) {
    try {
      synchronouslyWrite(mLeft, pair.r1());
      synchronouslyWrite(mRight, pair.r2());
    } catch (IOException e) {
      throw new RuntimeIOException(e);
    }
  }

  private void synchronouslyWrite(FastqWriter w, FastqSequence read) throws IOException {
    w.write(read.getName(), read.getBases(), read.getQualities(), read.getBases().length);
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
