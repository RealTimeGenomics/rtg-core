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

import com.rtg.calibrate.CalibrationStats;
import com.rtg.calibrate.Calibrator;
import com.rtg.calibrate.CovariateEnum;
import com.rtg.calibrate.StatsProcessor;
import com.rtg.reader.Arm;


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
    for (int i = 0; i < baseQualSize; ++i) {
      mCurve[i] = CalibratedMachineErrorParams.countsToEmpiricalQuality(proc.mMismatches[i], proc.mTotals[i], globalErrorRate);
//       System.out.println(i + ": miss=" + proc.mMismatches[i] + " total=" + proc.mTotals[i] + " qual=" + mCurve[i]);
    }
  }

  @Override
  public int getScaledPhred(byte quality, int readPosition, Arm arm) {
    final int qualIndex = quality & 0xFF;
    if (qualIndex >= mCurve.length) {
      return mCurve[mCurve.length - 1];
    }
    return mCurve[qualIndex];
  }

}

