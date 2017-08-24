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
}
