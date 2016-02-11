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

import java.util.Arrays;

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
    private final double[] mMeanSquareError;
    private final long[] mTerms;
    private final int mQualityCovariateIndex;

    private final int mReadPosIndex;

    public ReadPositionStatsProcessor(int readPosIndex, int readPosSize, int qualityCovariateIndex) {
      mMismatches = new long[readPosSize];
      mTotals = new long[readPosSize];
      mMeanSquareError = new double[readPosSize];
      mTerms = new long[readPosSize];
      mReadPosIndex = readPosIndex;
      mQualityCovariateIndex = qualityCovariateIndex;
    }

    @Override
    public void process(int[] covariateValues, CalibrationStats stats) {
      if (stats != null) {
        final int readPosValue = covariateValues[mReadPosIndex];
        final long different = stats.getDifferent();
        final long same = stats.getEqual();
        final long total = different + same;
        mMismatches[readPosValue] += different;
        mTotals[readPosValue] += total;
        if (mQualityCovariateIndex >= 0) {
          final int qual = stats.getCovariateValue(mQualityCovariateIndex);
          final int empirical = CalibratedMachineErrorParams.countsToEmpiricalQuality(different, total, 0); // todo purpose of global error rate
          final int error = qual - empirical;
          mMeanSquareError[readPosValue] += total * error * error; // todo is this total the right one to use?
          mTerms[readPosValue] += total;
        }
      }
    }
  }

  private final int[] mCurve;
  private final double[] mMeanSquareError;

  public ReadPositionPhredScaler(CovariateEnum type, Calibrator cal, Calibrator.QuerySpec query) {
    final int readPosIndex = cal.getCovariateIndex(type);
    final int readPosSize = cal.getCovariate(readPosIndex).size();
    final ReadPositionStatsProcessor proc = new ReadPositionStatsProcessor(readPosIndex, readPosSize, cal.getCovariateIndex(CovariateEnum.BASEQUALITY));
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
    mMeanSquareError = Arrays.copyOf(proc.mMeanSquareError, proc.mMeanSquareError.length); // XXX remove need for this
    for (int k = 0; k < readPosSize; k++) {
      mMeanSquareError[k] /= proc.mTerms[k];
    }
  }

  @Override
  public int getPhred(byte quality, int readPosition) {
    return mCurve[readPosition];
  }

  public double getMSE(final int readPosition) {
    return mMeanSquareError[readPosition];
  }

}
