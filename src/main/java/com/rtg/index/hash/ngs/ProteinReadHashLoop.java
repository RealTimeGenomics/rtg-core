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
