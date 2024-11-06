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
    Diagnostic.userLog("Calculated " + mDiscardTarget + "% volume hash frequency threshold: " + ret);
    final int preLimit = ret;
    if (mMaxThreshold > 0 && ret > mMaxThreshold) {
      ret = mMaxThreshold;
    }
    if (mMinThreshold > 0 && ret < mMinThreshold) {
      ret = mMinThreshold;
    }
    if (ret != preLimit) {
      Diagnostic.userLog("Limiting hash frequency threshold to bounds: [" + mMinThreshold + "," + mMaxThreshold + "] now " + ret);
    }
    return ret - 1;
  }

  @Override
  public String toString() {
    return "ProportionalRepeatFrequency(" + mDiscardTarget + "%, " + mMinThreshold + ", " + mMaxThreshold + ")";
  }
}
