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
