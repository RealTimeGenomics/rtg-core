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

import com.rtg.calibrate.CalibrationStats;
import com.rtg.calibrate.Calibrator;
import com.rtg.calibrate.CovariateEnum;
import com.rtg.calibrate.StatsProcessor;
import com.rtg.ngs.Arm;


/**
 * Implementation that uses calibration files.
 */
class BaseQualityMachineCyclePhredScaler implements PhredScaler {

  private static class BaseQualityReadPositionStatsProcessor implements StatsProcessor {
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
        mMismatches[baseQualityValue][readPosValue] += stats.getDifferent();
        mTotals[baseQualityValue][readPosValue] += stats.getDifferent() + stats.getEqual();
      }
    }
  }

  private final int[][] mCurve;

  BaseQualityMachineCyclePhredScaler(Calibrator cal, Calibrator.QuerySpec query) {
    final int baseQualIndex = cal.getCovariateIndex(CovariateEnum.BASEQUALITY);
    final int readPosIndex = cal.getCovariateIndex(CovariateEnum.MACHINECYCLE);
    final int baseQualSize = cal.getCovariate(baseQualIndex).size();
    final int readPosSize = cal.getCovariate(readPosIndex).size();
    final BaseQualityReadPositionStatsProcessor proc = new BaseQualityReadPositionStatsProcessor(baseQualIndex, baseQualSize, readPosIndex, readPosSize);
    cal.processStats(proc, query);
    long totalMismatches = 0;
    for (final long[] v : proc.mMismatches) {
      for (final long vv : v) {
        totalMismatches += vv;
      }
    }
    long totalEverything = 0;
    for (final long[] v : proc.mTotals) {
      for (final long vv : v) {
        totalEverything += vv;
      }
    }

    mCurve = new int[baseQualSize][readPosSize];
    final double globalErrorRate = (double) totalMismatches / (double) totalEverything;
    for (int i = 0; i < baseQualSize; i++) {
      for (int j = 0; j < readPosSize; j++) {
        mCurve[i][j] = CalibratedMachineErrorParams.countsToEmpiricalQuality(proc.mMismatches[i][j], proc.mTotals[i][j], globalErrorRate);
      }
    }
  }

  @Override
  public int getScaledPhred(byte quality, int readPosition, Arm arm) {
    final int qualIndex = quality & 0xFF;
    if (qualIndex >= mCurve.length) {
      return mCurve[mCurve.length - 1][readPosition];
    }
    return mCurve[qualIndex][readPosition];
  }

}
