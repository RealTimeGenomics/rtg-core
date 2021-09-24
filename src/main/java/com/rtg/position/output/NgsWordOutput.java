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
