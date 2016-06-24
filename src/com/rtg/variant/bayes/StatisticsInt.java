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
      mAmbiguous++;
    }

    if (evidence.isReadPaired()) {
      if (evidence.isMated()) {
        mMatedCount++;
      } else {
        mUnmatedCount++;
      }
    }

    mCounts.increment(evidence, bestHyp, errorIncrement);
    mTotalCoverage++;
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
    sb.append(mCounts.toString());
    return sb.toString();
  }

}
