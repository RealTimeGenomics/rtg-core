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


/**
 * Implementation that uses calibration files.
 */
class BaseQualityPhredScaler implements PhredScaler {

  private static class BaseQualityStatsProcessor implements StatsProcessor {

    private final long[] mMismatches;
    private final long[] mTotals;

    private final int mBaseQualIndex;

    BaseQualityStatsProcessor(int baseQualIndex, int baseQualSize) {
      mMismatches = new long[baseQualSize];
      mTotals = new long[baseQualSize];
      mBaseQualIndex = baseQualIndex;
    }
    @Override
    public void process(int[] covariateValues, CalibrationStats stats) {
      if (stats != null) {
        final int index = covariateValues[mBaseQualIndex];
        mMismatches[index] += stats.getDifferent();
        mTotals[index] += stats.getDifferent() + stats.getEqual();
      }
    }
  }

  private final int[] mCurve;

  BaseQualityPhredScaler(Calibrator cal, Calibrator.QuerySpec query) {
    final int baseQualIndex = cal.getCovariateIndex(CovariateEnum.BASEQUALITY);
    final int baseQualSize = cal.getCovariate(baseQualIndex).size();
    final BaseQualityStatsProcessor proc = new BaseQualityStatsProcessor(baseQualIndex, baseQualSize);
    cal.processStats(proc, query);
    long totalMismatches = 0;
    for (final long v : proc.mMismatches) {
      totalMismatches += v;
    }
    long totalEverything = 0;
    for (final long v : proc.mTotals) {
      totalEverything += v;
    }
    mCurve = new int[baseQualSize];
    // TODO we should use a laplace style estimator here
    final double globalErrorRate = (double) totalMismatches / (double) totalEverything;
//     System.out.println("global miss=" + totalMismatches + " total=" + totalEverything + " errRate=" + globalErrorRate);
    for (int i = 0; i < baseQualSize; i++) {
      mCurve[i] = CalibratedMachineErrorParams.countsToEmpiricalQuality(proc.mMismatches[i], proc.mTotals[i], globalErrorRate);
//       System.out.println(i + ": miss=" + proc.mMismatches[i] + " total=" + proc.mTotals[i] + " qual=" + mCurve[i]);
    }
  }

  @Override
  public int getPhred(byte quality, int readPosition) {
    final int qualIndex = quality & 0xFF;
    if (qualIndex >= mCurve.length) {
      return mCurve[mCurve.length - 1];
    }
    return mCurve[qualIndex];
  }

}

