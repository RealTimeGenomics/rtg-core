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
package com.rtg.index;

/**
 * Use an absolute threshold on frequency of hash occurrence
 */
public class FixedRepeatFrequencyFilterMethod implements IndexFilterMethod {

  /**
   * Hashes with more than this number of hits are assumed to be repeats and are
   * not reported during searches. This can give performance improvements by
   * eliminating the construction of index hit objects that only contain these
   * high frequency chunks.
   */
  protected final int mThreshold;

  /**
   * @param threshold maximum repeat frequency threshold.
   */
  public FixedRepeatFrequencyFilterMethod(final int threshold) {
    if (threshold < 1) {
      throw new RuntimeException("Threshold must be positive:" + threshold);
    }
    mThreshold = threshold;
  }

  @Override
  public IndexFilterMethod threadClone() {
    return new FixedRepeatFrequencyFilterMethod(mThreshold);
  }

  @Override
  public void initialize(Index index) {
    // Nothing to do
  }

  @Override
  public boolean keepHash(long hash, long numHits) {
    return numHits <= mThreshold;
  }

  @Override
  public String toString() {
    return "FixedRepeatFrequency(" + mThreshold + ")";
  }
}
