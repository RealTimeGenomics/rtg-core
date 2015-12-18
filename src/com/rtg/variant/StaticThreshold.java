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
package com.rtg.variant;

/**
 * threshold object that will always return the same value
 *         Date: 25/10/11
 *         Time: 9:52 AM
 */
public class StaticThreshold implements CoverageThreshold {
  private final int mThreshold;

  private final int mTotalThreshold;

  /**
   * Construct a threshold object that will always return <code>threshold</code>
   * @param threshold the threshold to return when asked
   */
  public StaticThreshold(int threshold) {
    this(threshold, threshold);
  }

  /**
   * @param singleThreshold threshold to use for {@link #thresholdSingle(String)}
   * @param totalThreshold threshold to use for {@link #thresholdTotal(String)}
   */
  public StaticThreshold(int singleThreshold, int totalThreshold) {
    mThreshold = singleThreshold;
    mTotalThreshold = totalThreshold;
  }

  @Override
  public int thresholdSingle(String sequenceName) {
    return mThreshold;
  }

  @Override
  public int thresholdTotal(String sequenceName) {
    return mTotalThreshold;
  }

  @Override
  public String toString() {
    return "" + mThreshold + ":" + mTotalThreshold;
  }

}
