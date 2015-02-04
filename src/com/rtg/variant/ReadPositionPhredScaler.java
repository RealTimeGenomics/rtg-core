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
 * Implementation that differentiates on read position only (useful
 * for overall analysis of average quality per position in the read).
 */
class ReadPositionPhredScaler implements PhredScaler {

  private static class ReadPositionStatsProcessor implements StatsProcessor {
    private final long[] mMismatches;
    private final long[] mTotals;

    private final int mReadPosIndex;

    public ReadPositionStatsProcessor(int readPosIndex, int readPosSize) {
      mMismatches = new long[readPosSize];
      mTotals = new long[readPosSize];
      mReadPosIndex = readPosIndex;
    }

    @Override
    public void process(int[] covariateValues, CalibrationStats stats) {
      if (stats != null) {
        final int readPosValue = covariateValues[mReadPosIndex];
        mMismatches[readPosValue] += stats.getDifferent();
        mTotals[readPosValue] += stats.getDifferent() + stats.getEqual();
      }
    }
  }

  private final int[] mCurve;

  public ReadPositionPhredScaler(Calibrator cal, Calibrator.QuerySpec query) {
    final int readPosIndex = cal.getCovariateIndex(CovariateEnum.READPOSITION);
    final int readPosSize = cal.getCovariate(readPosIndex).size();
    final ReadPositionStatsProcessor proc = new ReadPositionStatsProcessor(readPosIndex, readPosSize);
    cal.processStats(proc, query);
    long totalMismatches = 0;
    for (final long v : proc.mMismatches) {
      totalMismatches += v;
    }
    long totalEverything = 0;
    for (final long v : proc.mTotals) {
      totalEverything += v;
    }

    mCurve = new int[readPosSize];
    final double globalErrorRate = (double) totalMismatches / (double) totalEverything;
    for (int i = 0; i < readPosSize; i++) {
      mCurve[i] = CalibratedMachineErrorParams.countsToEmpiricalQuality(proc.mMismatches[i], proc.mTotals[i], globalErrorRate);
    }
  }

  @Override
  public int getPhred(byte quality, int readPosition) {
    return mCurve[readPosition];
  }

}
