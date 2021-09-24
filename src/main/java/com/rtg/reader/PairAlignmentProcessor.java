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

/**
 */
class PairAlignmentProcessor implements Runnable {
  private final PairAlignmentStats mStats;
  private final BatchReorderingWriter<FastqPair> mWriter;
  private final Batch<FastqPair> mBatch;
  private final PairAligner mAligner;

  PairAlignmentProcessor(PairAlignmentStats stats, BatchReorderingWriter<FastqPair> writer, Batch<FastqPair> batch, PairAligner aligner) {
    mStats = stats;
    mWriter = writer;
    mBatch = batch;
    mAligner = aligner;
  }


  @Override
  public void run() {
    mAligner.processReads(mBatch.getItems());
    mWriter.writeBatch(mBatch);
    synchronized (mStats) {
      mStats.accumulate(mAligner.getStats());
    }
  }
}
