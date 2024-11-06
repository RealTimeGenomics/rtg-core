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
