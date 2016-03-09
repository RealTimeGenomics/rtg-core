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
  /** do you want to interpolate from qualities out of the observer ranges per position? */
  private static final long QUALITY_MIN_EVIDENCE = 1000;

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
        if (proc.mTotals[i][j] < QUALITY_MIN_EVIDENCE) {
          // Save it for interpolation
          mCurve[i][j] = -1;
        } else {
          mCurve[i][j] = CalibratedMachineErrorParams.countsToEmpiricalQuality(proc.mMismatches[i][j], proc.mTotals[i][j], globalErrorRate);
        }
      }
    }

    //Interpolation phase
    // First work out read position independent quality levels to be used as defaults
    final int[] qualityMismatches = new int[mCurve.length];
    final int[] qualityTotals = new int[mCurve.length];
    for (int i = 0; i < qualityMismatches.length; i++) {
      for (int j = 0; j < proc.mTotals[i].length; j++) {
        qualityMismatches[i] += proc.mMismatches[i][j];
        qualityTotals[i] += proc.mTotals[i][j];
      }
    }

    final int[] qualityDefaults = new int[mCurve.length];
    for (int i = 0; i < qualityDefaults.length; i++) {
      if (qualityTotals[i] < QUALITY_MIN_EVIDENCE) {
        qualityDefaults[i] = -1;
      } else {
        qualityDefaults[i] = CalibratedMachineErrorParams.countsToEmpiricalQuality(qualityMismatches[i], qualityTotals[i], globalErrorRate);
      }
    }
    // Some of the quality default values will be missing so need to interpolate those as well;
    // if min/max qual are missing then interpolate between claimed and the first/last observed qual
    if (qualityDefaults[0] == -1) {
      qualityDefaults[0] = 0;
    }

    if (qualityDefaults[qualityDefaults.length - 1] == -1) {
      qualityDefaults[qualityDefaults.length - 1] = qualityDefaults.length - 1;
    }

    // Now we have defaults for the ends of curves we can interpolate within the read position curve set.
    // At each read position we'll draw a straight line between present quality values we have a rate for.
    // Use the global quality default array for ends...
    for (int j = 0; j < readPosSize; j++) {
      final Interpolate2dArrayColumn interpolate2dArrayColumn = new Interpolate2dArrayColumn(j, mCurve);
      interpolate2dArrayColumn.fill(qualityDefaults);
      interpolate2dArrayColumn.process();
    }
  }

  @Override
  public int getScaledPhred(byte quality, int readPosition, Arm arm) {
    final int qualIndex = quality & 0xFF;
    if (qualIndex >= mCurve.length) {
      return mCurve[mCurve.length - 1][Math.min(readPosition, mCurve[mCurve.length - 1].length - 1)];
    } else if (readPosition >= mCurve[qualIndex].length) {
      return mCurve[qualIndex][mCurve[qualIndex].length - 1];
    }
    return mCurve[qualIndex][readPosition];

  }

}
