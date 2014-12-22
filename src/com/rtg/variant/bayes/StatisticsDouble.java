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
package com.rtg.variant.bayes;

import com.rtg.util.MathUtils;
import com.rtg.util.StringUtils;

/**
 */
public class StatisticsDouble extends Statistics<AlleleStatisticsDouble> {

  private double mAmbiguous;

  private double mTotalCoverage = 0;
  private double mTotalError = 0;

  private double mMatedCount;
  private double mUnmatedCount;

  /**
   * @param description about which statistics are being collected.
   */
  public StatisticsDouble(Description description) {
    super(description, new AlleleStatisticsDouble(description));
  }

  @Override
  protected void incrementBest(final EvidenceInterface distribution, int bestHyp) {

    final double coverageIncrement = coverageIncrement(distribution);
    final double errorIncrement = errorIncrement(distribution);

    if (distribution.mapError() >= Model.AMBIGUITY_THRESHOLD) {
      mAmbiguous += coverageIncrement;
    }
    if (distribution.isReadPaired()) {
      if (distribution.isMated()) {
        mMatedCount += coverageIncrement;
      } else {
        mUnmatedCount += coverageIncrement;
      }
    }

    mCounts.increment(distribution, bestHyp, errorIncrement, coverageIncrement);
    mTotalCoverage += coverageIncrement;
    mTotalError += errorIncrement;
  }

  /** Return the amount that error statistics should be incremented by for this piece of evidence */
  protected double errorIncrement(EvidenceInterface distribution) {
    final double r = distribution.mapError();
    final double q = distribution.error();
    return r + (1.0 - r) * q;
  }

  /** Return the amount that coverage statistics should be incremented by for this piece of evidence */
  protected double coverageIncrement(EvidenceInterface distribution) {
    return 1.0;
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
    return (int) MathUtils.round(mTotalCoverage);
  }


  /**
   * @return the number of ambiguous reads.
   */
  @Override
  public final int ambiguousCount() {
    return (int) MathUtils.round(mAmbiguous);
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
    sb.append("coverage=").append(Double.toString(mTotalCoverage));
    final String correction = String.format("%1$04.3f", mTotalError);
    sb.append(" correction=").append(correction);
    sb.append(StringUtils.LS);
    sb.append(mCounts.toString());
    return sb.toString();
  }

}
