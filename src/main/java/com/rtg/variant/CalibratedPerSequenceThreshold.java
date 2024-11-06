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

import java.util.Collection;
import java.util.HashMap;
import java.util.Map;

import com.rtg.calibrate.CalibratedPerSequenceExpectedCoverage;
import com.rtg.util.MathUtils;
import com.rtg.util.diagnostic.NoTalkbackSlimException;

/**
 * Threshold that returns a multiple of the average coverage for the current sequence.
 */
public class CalibratedPerSequenceThreshold implements CoverageThreshold {

  // Thresholds per sequence, selected using maximum coverage across samples
  final Map <String, Integer> mMaxThresholds;
  // Thresholds per sequence, selected using total coverage across samples
  final Map <String, Integer> mSumThresholds;

  final double mMultiple;
  final ThresholdFunction mThresholdFunction;

  /**
   * Contains some alternative functions for obtaining a coverage threshold given
   * a coverage and a user-variable parameter.
   */
  public enum ThresholdFunction {

    /** Threshold is coverage * dial */
    SIMPLE_MULT {
      @Override
      public int threshold(double coverage, double dial) {
        return (int) MathUtils.round(coverage * dial);
      }
      @Override
      public String functionString(double dial) {
        return "avgSeqCov*" + dial;
      }
    },

    /**
     * An alternative average coverage multiplier based coverage threshold. John suggests
     * that this method scales better with varying average coverage than a simple multiplier.
     *
     * threshold = average_coverage + dial * sqrt(average_coverage)
     */
    SQRT_MULT {
      @Override
      public int threshold(double coverage, double dial) {
        return (int) MathUtils.round(coverage + Math.sqrt(coverage) * dial);
      }
      @Override
      public String functionString(double dial) {
        return "avgSeqCov+(sqrt(avgSeqCov)*" + dial + ")";
      }
    };

    /**
     * Return the coverage threshold given the average sequence coverage and a
     * user specified parameter.
     * @param coverage the average sequence coverage
     * @param dial the user-specified parameter
     * @return the coverage threshold
     */
    public abstract int threshold(double coverage, double dial);

    /**
     * Get the string of the function applied to get the coverage threshold.
     * @param dial the user-specified parameter
     * @return the string representing the mathematical function applied to get threshold.
     */
    public abstract String functionString(double dial);
  }

  /**
   *
   * @param coverages the mapping calibration stats to use when computing coverage
   * @param multiple factor to multiply coverage by when computing threshold
   * @param tfunc the <code>ThresholdFunction</code> to use
   */
  public CalibratedPerSequenceThreshold(CalibratedPerSequenceExpectedCoverage coverages, double multiple, ThresholdFunction tfunc) {
    mMultiple = multiple;
    mThresholdFunction = tfunc;
    final Collection<String> sequences = coverages.sequences();
    mMaxThresholds = new HashMap<>(sequences.size());
    mSumThresholds = new HashMap<>(sequences.size());
    for (final String seqName : sequences) {
      final Double sum = coverages.expectedTotalCoverage(seqName);
      mSumThresholds.put(seqName, mThresholdFunction.threshold(Math.max(1, sum), mMultiple));

      double currentMax = 0;
      for (String sample : coverages.samples()) {
        final double sampleCov = coverages.expectedCoverage(seqName, sample);
        if (sampleCov > currentMax) {
          currentMax = sampleCov;
        }
      }
      mMaxThresholds.put(seqName, mThresholdFunction.threshold(Math.max(1, currentMax), mMultiple));
    }
  }

  @Override
  public int thresholdSingle(String sequenceName) {
    if (!mMaxThresholds.containsKey(sequenceName)) {
      throw new NoTalkbackSlimException("Unknown sequence: " + sequenceName);
    }
    return mMaxThresholds.get(sequenceName);
  }

  @Override
  public int thresholdTotal(String sequenceName) {
    if (!mSumThresholds.containsKey(sequenceName)) {
      throw new NoTalkbackSlimException("Unknown sequence: " + sequenceName);
    }
    return mSumThresholds.get(sequenceName);
  }

  @Override
  public String toString() {
    return mThresholdFunction.functionString(mMultiple);
  }

}
