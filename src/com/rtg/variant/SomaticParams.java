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

import com.rtg.util.StringUtils;
import com.rtg.util.intervals.ReferenceRanges;

/**
 * @author kurt
 */
public class SomaticParams {

  private final double mSomaticRate;
  private final boolean mIncludeGermlineVariants;
  private final boolean mIncludeGainOfReference;
  private final ReferenceRanges<Double> mSiteSpecificSomaticPriors;
  private final double mLohPrior;

  /**
   * @param builder the builder object.
   */
  SomaticParams(SomaticParamsBuilder builder) {
    mSomaticRate = builder.mSomaticRate;
    mIncludeGermlineVariants = builder.mIncludeGermlineVariants;
    mIncludeGainOfReference = builder.mIncludeGainOfReference;
    mSiteSpecificSomaticPriors = builder.mSiteSpecificSomaticPriors;
    mLohPrior = builder.mLohPrior;
  }

  /**
   * @return expected mutation rate in the cancer cells with respect to somatic cells.
   */
  public double somaticRate() {
    return mSomaticRate;
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
   * @return prior probability that a region in a cancer genome has lost heterozygosity.
   */
  public double lohPrior() {
    return mLohPrior;
  }

  @Override
  public String toString() {
    final StringBuilder sb = new StringBuilder(super.toString());
    sb.append("VariantParams")
      .append(" somaticRate=").append(somaticRate())
      .append(" includeGermlineVariants=").append(includeGermlineVariants())
      .append(" includeGainOfReference=").append(includeGainOfReference())
      .append(" lohPrior=").append(lohPrior()).append(StringUtils.LS);
    return sb.toString();
  }
}
