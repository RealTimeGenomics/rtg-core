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
