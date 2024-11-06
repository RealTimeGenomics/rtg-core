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

package com.rtg.variant;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.calibrate.CalibrationStats;
import com.rtg.calibrate.StatsProcessor;

/**
 *
 */
@TestClass("com.rtg.variant.BaseQualityMachineCyclePhredScalerTest")
class BaseQualityReadPositionStatsProcessor implements StatsProcessor {
  private final long[][] mMismatches;
  private final long[][] mTotals;

  private final int mBaseQualIndex;
  private final int mReadPosIndex;

  BaseQualityReadPositionStatsProcessor(int baseQualIndex, int baseQualSize, int readPosIndex, int readPosSize) {
    mMismatches = new long[baseQualSize][readPosSize];
    mTotals = new long[baseQualSize][readPosSize];
    mBaseQualIndex = baseQualIndex;
    mReadPosIndex = readPosIndex;
  }

  @Override
  public void process(int[] covariateValues, CalibrationStats stats) {
    if (stats != null) {
      final int baseQualityValue = covariateValues[mBaseQualIndex];
      final int readPosValue = covariateValues[mReadPosIndex];
      getMismatches()[baseQualityValue][readPosValue] += stats.getDifferent();
      getTotals()[baseQualityValue][readPosValue] += stats.getDifferent() + stats.getEqual();
    }
  }

  public long[][] getMismatches() {
    return mMismatches;
  }

  public long[][] getTotals() {
    return mTotals;
  }
}
