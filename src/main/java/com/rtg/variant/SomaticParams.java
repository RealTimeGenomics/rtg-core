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

import com.rtg.util.intervals.ReferenceRanges;

/**
 */
public class SomaticParams {

  private final double mSomaticRate;
  private final boolean mIncludeGermlineVariants;
  private final boolean mIncludeGainOfReference;
  private final ReferenceRanges<Double> mSiteSpecificSomaticPriors;
  private final double mLohPrior;
  private final boolean mSomaticAlleleBalance;
  private final int mContaminationBasis;

  /**
   * @param builder the builder object.
   */
  SomaticParams(SomaticParamsBuilder builder) {
    mSomaticRate = builder.mSomaticRate;
    mSomaticAlleleBalance = builder.mSomaticAlleleBalance;
    mIncludeGermlineVariants = builder.mIncludeGermlineVariants;
    mIncludeGainOfReference = builder.mIncludeGainOfReference;
    mSiteSpecificSomaticPriors = builder.mSiteSpecificSomaticPriors;
    mLohPrior = builder.mLohPrior;
    mContaminationBasis = builder.mContaminationBasis;
  }

  /**
   * @return expected mutation rate in the cancer cells with respect to somatic cells.
   */
  public double somaticRate() {
    return mSomaticRate;
  }

  /**
   * @return true if somatic allele balance computation should be employed
   */
  public boolean somaticAlleleBalance() {
    return mSomaticAlleleBalance;
  }

  /**
   * @return true iff germline variants should be included in VCF output during somatic calling.
   */
  public boolean includeGermlineVariants() {
    return mIncludeGermlineVariants;
  }

  /**
   * @return true iff gain of reference somatic calls should be output.
   */
  public boolean includeGainOfReference() {
    return mIncludeGainOfReference;
  }

  /**
   * @return site specific somatic priors.
   */
  public ReferenceRanges<Double> siteSpecificSomaticPriors() {
    return mSiteSpecificSomaticPriors;
  }

  /**
   * @return number of examples for contamination estimation
   */
  public int contaminationBasis() {
    return mContaminationBasis;
  }

  /**
   * @return prior probability that a region in a cancer genome has lost heterozygosity.
   */
  public double lohPrior() {
    return mLohPrior;
  }

  @Override
  public String toString() {
    final StringBuilder sb = new StringBuilder();
    sb.append("SomaticParams")
      .append(" somatic_rate=").append(somaticRate())
      .append(" somatic_allele_balance=").append(somaticAlleleBalance())
      .append(" somatic_ssp=").append(siteSpecificSomaticPriors() == null ? "none" : "supplied")
      .append(" include_germline_variants=").append(includeGermlineVariants())
      .append(" include_gain_of_reference=").append(includeGainOfReference())
      .append(" loh_prior=").append(lohPrior())
      .append(" contamination_basis=").append(contaminationBasis());
    return sb.toString();
  }
}
