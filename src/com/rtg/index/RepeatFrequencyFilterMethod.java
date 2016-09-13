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

import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.OneShotTimer;

/**
 *
 */
public class RepeatFrequencyFilterMethod implements IndexFilterMethod {

  protected final boolean mUseProportionalThreshold;

  /**
   * Hashes with more than this number of hits are assumed to be repeats and are
   * not reported during searches. This can give performance improvements by
   * eliminating the construction of index hit objects that only contain these
   * high frequency chunks.
   */
  protected final int mThreshold;
  /**
   * When using proportional threshold don't exceed this value as cutoff.
   */
  protected final int mMaxThreshold;
  /**
   * When using proportional threshold don't dont go below this value.
   */
  protected final int mMinThreshold;

  private int mCutoff;

  /**
   * @param threshold maximum repeat frequency threshold - default
   *        <code>Integer.MAX_VALUE</code> if null.
   * @param proportionalThreshold Whether the frequency threshold should be calculated from index data rather than as a parameter.
   * @param maxThreshold when using proportional threshold don't exceed this repeat frequency
   * @param minThreshold when using proportional threshold don't go below this repeat frequency
   */
  public RepeatFrequencyFilterMethod(final Integer threshold, boolean proportionalThreshold, int maxThreshold, int minThreshold) {
    if (threshold == null) {
      mThreshold = Integer.MAX_VALUE;
    } else {
      if (!proportionalThreshold && threshold < 1) {
        throw new RuntimeException("Threshold must be positive:" + threshold);
      }
      mThreshold = threshold;
    }
    mMaxThreshold = maxThreshold;
    mMinThreshold = minThreshold;

    mUseProportionalThreshold = proportionalThreshold;
  }

  @Override
  public IndexFilterMethod threadClone() {
    return new RepeatFrequencyFilterMethod(mThreshold, mUseProportionalThreshold, mMaxThreshold, mMinThreshold);
  }

  @Override
  public void initialize(Index index) {
    if (mUseProportionalThreshold) {
      final OneShotTimer over = new OneShotTimer("Index_Frequency");
      internalInitializeProportional(index.getSparseFrequencyHistogram(), index.getInitialHashes());
      over.stopLog();
    } else {
      mCutoff = mThreshold;
    }
  }

  void internalInitializeProportional(SparseFrequencyHistogram histo, long initialHashes) {
    mCutoff = determineCutoffFrequency(histo, initialHashes, mThreshold);
  }

  @Override
  public boolean keepHash(long hash, long numHits) {
    return numHits <= mCutoff;
  }

  protected int determineCutoffFrequency(final SparseFrequencyHistogram freqHist, final long initialHashes, final int percent) {
    final long cutoff = initialHashes * percent / 100;
    long val = 0;
    int ret = Integer.MAX_VALUE;
    for (int i = freqHist.length() - 1; i >= 0; i--) {
      final int freq = freqHist.getFrequency(i);
      val += freqHist.getCount(i) * freq;
      if (val > cutoff) {
        break;
      }
      ret = freq;
    }
    Diagnostic.userLog("Calculated hash frequency threshold: " + ret);
    if (mUseProportionalThreshold && (ret > mMaxThreshold || ret < mMinThreshold)) {
      ret = Math.max(mMinThreshold, Math.min(ret, mMaxThreshold));
      Diagnostic.userLog("Limiting hash frequency to bounds: [" + mMinThreshold + "," + mMaxThreshold + "] now " + ret);
    }
    return ret - 1;
  }
}
