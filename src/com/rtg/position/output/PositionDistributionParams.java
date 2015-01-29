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
package com.rtg.position.output;

import com.rtg.util.ObjectParams;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.Integrity;

/**
 * Specifies output parameters for Position.
 */
public class PositionDistributionParams extends ObjectParams implements Integrity {

  private final double mSubsProb;

  private final double mIndelOpenProb;

  private final double mIndelExtendProb;

  private final Integer mMaxGap;

  private final int mMaxIndel;

  /**
   * Use this version with Blosum62 scoring.
   * @param maxGap maximum gap permitted in map and region formats.
   */
  public PositionDistributionParams(final Integer maxGap) {
    this(Double.NaN,  Double.NaN, maxGap, 0);
  }


  /**
   * @param subsProb probability of a substitution (mismatch).
   * @param indelOpenProb probability that a new indel is to be opened.
   * @param maxGap maximum gap permitted in map and region formats.
   * @param maxIndel max indel
   */
  public PositionDistributionParams(final double subsProb, final double indelOpenProb, final Integer maxGap, final int maxIndel) {
    mSubsProb = subsProb;
    mIndelOpenProb = indelOpenProb;
    mIndelExtendProb = indelOpenProb;
    mMaxGap = maxGap;
    mMaxIndel = maxIndel;
    mObjects = new Object[] {mSubsProb, mIndelOpenProb, mIndelExtendProb, mMaxGap};
  }

  /**
   * Get <code>subsProb</code>.
   * @return Returns the substitution probability.
   */
  public double subsProb() {
    return mSubsProb;
  }

  /**
   * Get <code>indelOpenProb</code>.
   * @return Returns the indel open probability.
   */
  public double indelOpenProb() {
    return mIndelOpenProb;
  }

  /**
   * Get <code>indelExtendProb</code>.
   * @return Returns the indel extend probability.
   */
  public double indelExtendProb() {
    return mIndelExtendProb;
  }

  /**
   * Get the maximum allowed gap.
   * @return the maximum allowed gap.
   */
  public Integer maxGap() {
    return mMaxGap;
  }

  /**
   * When using long reads get max indel size
   * @return the value
   */
  public int maxIndel() {
    return mMaxIndel;
  }

  @Override
  public void close() {
    // do nothing
  }

  @Override
  public String toString() {
    final StringBuilder sb = new StringBuilder();
    if (!Double.isNaN(mSubsProb) || !Double.isNaN(mIndelOpenProb) || !Double.isNaN(mIndelExtendProb)) {
      sb.append("subs=").append(mSubsProb).append(" indelOpen=").append(mIndelOpenProb).append(" indelExtend=").append(mIndelExtendProb).append(" ");
    }
    sb.append("maxGap=").append(mMaxGap);
    return sb.toString();
  }

  @Override
  public boolean integrity() {
    if (Double.isNaN(mSubsProb)) {
      Exam.assertTrue(Double.isNaN(mIndelExtendProb));
      Exam.assertTrue(Double.isNaN(mIndelOpenProb));
    } else {
      Exam.assertTrue(GappedDistribution.prob(mIndelExtendProb));
      Exam.assertTrue(GappedDistribution.prob(mIndelOpenProb));
    }
    if (mMaxGap != null) {
      Exam.assertTrue("mMaxGap=" + mMaxGap, mMaxGap > 0);
    }
    Exam.assertTrue(mMaxIndel >= 0);
    return true;
  }

  @Override
  public boolean globalIntegrity() {
    integrity();
    return true;
  }
}
