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

import com.rtg.calibrate.Calibrator;
import com.rtg.calibrate.CovariateEnum;
import com.rtg.launcher.GlobalFlags;
import com.rtg.ngs.Arm;


/**
 * Implementation that uses calibration files.
 */
class CovariateIntersectionCycleQualityScaler implements PhredScaler {

  private final int[][] mCurve;

  CovariateIntersectionCycleQualityScaler(Calibrator cal, Calibrator.QuerySpec query) {
    final int baseQualIndex = cal.getCovariateIndex(CovariateEnum.BASEQUALITY);
    final int readPosIndex = cal.getCovariateIndex(CovariateEnum.MACHINECYCLE);
    final int baseQualSize = cal.getCovariate(baseQualIndex).size();
    final int readPosSize = cal.getCovariate(readPosIndex).size();
    final BaseQualityReadPositionStatsProcessor proc = new BaseQualityReadPositionStatsProcessor(baseQualIndex, baseQualSize, readPosIndex, readPosSize);
    cal.processStats(proc, query);
    long totalMismatches = 0;
    for (final long[] v : proc.getMismatches()) {
      for (final long vv : v) {
        totalMismatches += vv;
      }
    }
    long totalEverything = 0;
    for (final long[] v : proc.getTotals()) {
      for (final long vv : v) {
        totalEverything += vv;
      }
    }

    final double globalErrorRate = (double) totalMismatches / (double) totalEverything;
    final int globalPhred = CalibratedMachineErrorParams.countsToEmpiricalQuality(totalMismatches, totalEverything, globalErrorRate);

    final long[] qualityErrors = new long[baseQualSize];
    final long[] cycleErrors = new long[readPosSize];
    final long[] qualityTotals = new long[baseQualSize];
    final long[] cycleTotals = new long[readPosSize];
    for (int i = 0; i < baseQualSize; i++) {
      for (int j = 0; j < readPosSize; j++) {
        qualityTotals[i] += proc.getTotals()[i][j];
        qualityErrors[i] += proc.getMismatches()[i][j];
        cycleTotals[j] += proc.getTotals()[i][j];
        cycleErrors[j] += proc.getMismatches()[i][j];
      }
    }

    mCurve = new int[baseQualSize][readPosSize];
    final long qualityCalibrationMinEvidence = GlobalFlags.getIntegerValue(GlobalFlags.QUALITY_CALIBRATION_MIN_EVIDENCE);
    for (int i = 0; i < baseQualSize; i++) {
      for (int j = 0; j < readPosSize; j++) {
        if (proc.getMismatches()[i][j] < qualityCalibrationMinEvidence) {
          // Save it for interpolation
          mCurve[i][j] = -1;
        } else {
          mCurve[i][j] = CalibratedMachineErrorParams.countsToEmpiricalQuality(proc.getMismatches()[i][j], proc.getTotals()[i][j], globalErrorRate);
        }
      }
    }

    for (int i = 0; i < baseQualSize; i++) {
      for (int j = 0; j < readPosSize; j++) {
        if (mCurve[i][j] == -1 && qualityErrors[i] > qualityCalibrationMinEvidence && cycleErrors[j] > qualityCalibrationMinEvidence) {
          final int qualityAdjust = CalibratedMachineErrorParams.countsToEmpiricalQuality(qualityErrors[i], qualityTotals[i], globalErrorRate) - globalPhred;
          final int positionAdjust = CalibratedMachineErrorParams.countsToEmpiricalQuality(cycleErrors[j], cycleTotals[j], globalErrorRate) - globalPhred;
          // >>> 1 to convince findbugs this is overflow safe
          mCurve[i][j] = (globalPhred * 2 + qualityAdjust + positionAdjust) >>> 1;
        }
      }
    }

    //Interpolation phase
    // First work out read position independent quality levels to be used as defaults
    final int[] qualityDefaults = new int[mCurve.length];
    for (int i = 0; i < qualityDefaults.length; i++) {
      if (qualityTotals[i] < qualityCalibrationMinEvidence) {
        qualityDefaults[i] = -1;
      } else {
        qualityDefaults[i] = CalibratedMachineErrorParams.countsToEmpiricalQuality(qualityErrors[i], qualityTotals[i], globalErrorRate);
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
      final Interpolate2dArray interpolate2DArray = Interpolate2dArray.column(mCurve, j);
      interpolate2DArray.fill(qualityDefaults);
      interpolate2DArray.interpolate();
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
