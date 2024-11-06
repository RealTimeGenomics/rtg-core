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

package com.rtg.sam;

import java.io.IOException;
import java.util.Iterator;

import com.rtg.util.Populator;
import com.rtg.util.diagnostic.SpyTimer;
import com.rtg.util.diagnostic.Timer;

/**
 * Thread safe access to everything that can be used concurrently in {@link com.rtg.variant.bayes.multisample.MultisampleTask}.
 * @param <T> record type
 */
public final class CircularBufferMultifileSinglePassReaderWindowSync<T extends ReaderRecord<T> & MateInfo> extends CircularBufferMultifileSinglePassReaderWindow<T> {

  private static final Timer FLUSH_TIMER = new SpyTimer("flush");
  private static final Timer RECORDS_OVERLAP_TIMER = new SpyTimer("recordsOverlap");

  /**
   * Note: the supplied record iterator must be closed explicitly by the caller when finished using it as the
   * close method for the <code>CircularBufferMultifileSinglePassReaderWindow</code> will not close it.
   * @param recordIt iterator supplying records
   * @param pop the record populator
   * @param templateIndex sequence index for template
   * @param templateStart first position being retrieved
   * @param maxDepth maximum records per position
   */
  public CircularBufferMultifileSinglePassReaderWindowSync(RecordIterator<T> recordIt, Populator<T> pop, int templateIndex, int templateStart, int maxDepth) {
    super(recordIt, pop, templateIndex, templateStart, maxDepth);
  }

  @Override
  public synchronized void flush(int start, int end) throws IOException {
    //System.err.println("CirularBuffer flush enter");
    FLUSH_TIMER.start();
    try {
      super.flush(start, end);
    } finally {
      FLUSH_TIMER.stop();
    }
    //System.err.println("CirularBuffer flush exit");
  }

  @Override
  public synchronized Iterator<T> recordsOverlap(int start, int end) throws IOException {
    //System.err.println("CirularBuffer recordsOverlap enter");
    RECORDS_OVERLAP_TIMER.start();
    try {
      return super.recordsOverlap(start, end);
    } finally {
      RECORDS_OVERLAP_TIMER.stop();
      //System.err.println("CirularBuffer recordsOverlap exit");
    }
  }

  @Override
  public synchronized void advanceBuffer(int end) {
    super.advanceBuffer(end);
  }

  @Override
  public synchronized int flushedTo() {
    return super.flushedTo();
  }

  @Override
  public synchronized int finishedTo() {
    return super.finishedTo();
  }
}
