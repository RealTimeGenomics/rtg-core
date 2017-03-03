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
package com.rtg.alignment;

/**
 * Adapts a Unidirectional edit distance to an omni directional one, with the understanding that all calls are to be non
 * reverse complement.
 */
public class UnidirectionalAdaptor implements EditDistance {

  final UnidirectionalEditDistance mEd;

  /**
   * @param ed a unidirectional edit distance to adapt to an omnidirectional edit distance
   */
  public UnidirectionalAdaptor(UnidirectionalEditDistance ed) {
    mEd = ed;
  }

  @Override
  public int[] calculateEditDistance(byte[] read, int rlen, byte[] template, int zeroBasedStart, boolean rc, int maxScore, int maxShift, boolean cgLeft) {
    assert !rc;
    return mEd.calculateEditDistance(read, rlen, template, zeroBasedStart, maxScore, maxShift, cgLeft);
  }

  @Override
  public void logStats() {
    mEd.logStats();
  }
}
