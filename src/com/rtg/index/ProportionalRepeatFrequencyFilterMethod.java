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
 * Discards hashes occurring more frequently than a threshold determined to target a percentage of the hash volume.
 */
public class ProportionalRepeatFrequencyFilterMethod implements IndexFilterMethod {

  /**
   * Target percentage of hash volume to discard
   */
  protected final int mDiscardTarget;
  /**
   * When selecting the frequency threshold, don't exceed this value.
   */
  protected final int mMaxThreshold;
  /**
   * When selecting the frequency threshold, don't go below this value.
   */
  protected final int mMinThreshold;

  /**
   * Hashes with more than this number of hits are assumed to be repeats and are
   * not reported during searches. This can give performance improvements by
   * eliminating the construction of index hit objects that only contain these
   * high frequency chunks.
   */
  private int mCutoff;

  /**
   * @param discardTarget target percentage of hash volume to discard
   * @param maxThreshold when using proportional threshold don't exceed this repeat frequency, -1 to disable bounding
   * @param minThreshold when using proportional threshold don't go below this repeat frequency, -1 to disable bounding
   */
  public ProportionalRepeatFrequencyFilterMethod(final int discardTarget, int maxThreshold, int minThreshold) {
    mDiscardTarget = discardTarget;
    mMaxThreshold = maxThreshold;
    mMinThreshold = minThreshold;
  }

  @Override
  public IndexFilterMethod threadClone() {
    return new ProportionalRepeatFrequencyFilterMethod(mDiscardTarget, mMaxThreshold, mMinThreshold);
  }

  @Override
  public void initialize(Index index) {
    final OneShotTimer over = new OneShotTimer("Index_Frequency");
    internalInitializeProportional(index.getSparseFrequencyHistogram(), index.getInitialHashes());
    over.stopLog();
  }

  void internalInitializeProportional(SparseFrequencyHistogram histo, long initialHashes) {
    mCutoff = determineCutoffFrequency(histo, initialHashes, mDiscardTarget);
  }

  @Override
  public boolean keepHash(long hash, long numHits) {
    return numHits <= mCutoff;
  }

  private int determineCutoffFrequency(final SparseFrequencyHistogram freqHist, final long initialHashes, final int percent) {
    final long cutoff = initialHashes * percent / 100;
    long val = 0;
    int ret = Integer.MAX_VALUE;
    for (int i = freqHist.length() - 1; i >= 0; --i) {
      final int freq = freqHist.getFrequency(i);
      val += freqHist.getCount(i) * freq;
      if (val > cutoff) {
        break;
      }
      ret = freq;
    }
    Diagnostic.userLog("Calculated hash frequency threshold: " + ret);
    final int preLimit = ret;
    if (mMaxThreshold > 0 && ret > mMaxThreshold) {
      ret = mMaxThreshold;
    }
    if (mMinThreshold > 0 && ret < mMinThreshold) {
      ret = mMinThreshold;
    }
    if (ret != preLimit) {
      Diagnostic.userLog("Limiting hash frequency to bounds: [" + mMinThreshold + "," + mMaxThreshold + "] now " + ret);
    }
    return ret - 1;
  }
}
