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

package com.rtg.variant.bayes;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.util.StringUtils;

/**
 * Maintains various counts etc. that are used by the model logic during output.
 * It is expected that there will be classes of this that maintain other information that
 * is presented in human readable form.
 */
@TestClass("com.rtg.variant.bayes.snp.StatisticsSnpTest")
public class StatisticsInt extends Statistics<AlleleStatisticsInt> {

  private int mAmbiguous;

  private int mTotalCoverage;
  private double mTotalError;
  private int mMatedCount;
  private int mUnmatedCount;

  /**
   * @param description about which statistics are being collected.
   */
  public StatisticsInt(Description description) {
    super(description, new AlleleStatisticsInt(description));
  }

  @Override
  protected void incrementBest(final EvidenceInterface evidence, int bestHyp) {

    final double r = evidence.mapError();
    final double q = evidence.error();
    final double errorIncrement = r + (1.0 - r) * q;

    if (evidence.mapError() >= Model.AMBIGUITY_THRESHOLD) {
      ++mAmbiguous;
    }

    if (evidence.isReadPaired()) {
      if (evidence.isMated()) {
        ++mMatedCount;
      } else {
        ++mUnmatedCount;
      }
    }

    mCounts.increment(evidence, bestHyp, errorIncrement);
    ++mTotalCoverage;
    mTotalError += errorIncrement;
  }

  @Override
  protected double unmatedProbability() {
    final double total = mMatedCount + mUnmatedCount;
    if (total > 0) {
      return mUnmatedCount / total;
    }
    return 0.0;
  }

  /**
   * @return the effective coverage for the model
   */
  @Override
  public int coverage() {
    return mTotalCoverage;
  }

  /**
   * @return the number of ambiguous reads.
   */
  @Override
  public final int ambiguousCount() {
    return mAmbiguous;
  }

  /**
   * @return the total of the error corrections.
   */
  @Override
  public double totalError() {
    return mTotalError;
  }

  @Override
  public String toString() {
    // only used for debugging
    final StringBuilder sb = new StringBuilder();
    sb.append("Counts ");
    sb.append(StringUtils.LS);
    sb.append("coverage=").append(Integer.toString(mTotalCoverage));
    final String correction = String.format("%1$04.3f", mTotalError);
    sb.append(" correction=").append(correction);
    sb.append(StringUtils.LS);
    sb.append(mCounts);
    return sb.toString();
  }

}
