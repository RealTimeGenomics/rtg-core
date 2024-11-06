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

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.util.intervals.ReferenceRanges;

/**
 */
@TestClass("com.rtg.variant.SomaticParamsTest")
public class SomaticParamsBuilder {
  double mSomaticRate = 0.3;
  double mLohPrior = 0.0;
  boolean mIncludeGermlineVariants = false;
  boolean mIncludeGainOfReference = false;
  ReferenceRanges<Double> mSiteSpecificSomaticPriors = null;
  boolean mSomaticAlleleBalance = false;
  int mContaminationBasis = 1000;

  /**
   * Set the somatic rate
   * @param p somatic rate
   * @return this builder, so calls can be chained.
   */
  public SomaticParamsBuilder somaticRate(final double p) {
    if (p <= 0 || p >= 1) {
      throw new IllegalArgumentException();
    }
    mSomaticRate = p;
    return self();
  }

  /**
   * If set, incorporate expected somatic allelic fraction into scoring.
   * @param somaticAlleleBalance true iff somatic allele balance computation should be used.
   * @return this, for chaining
   */
  public SomaticParamsBuilder somaticAlleleBalance(boolean somaticAlleleBalance) {
    mSomaticAlleleBalance = somaticAlleleBalance;
    return self();
  }

  /**
   * If set, output germline variants in addition to somatic variants during somatic calling.
   * @param includeGermline true iff germline variants should be output
   * @return this, for chaining
   */
  public SomaticParamsBuilder includeGermlineVariants(boolean includeGermline) {
    mIncludeGermlineVariants = includeGermline;
    return self();
  }

  /**
   * If set, output should include somatic calls where there has been a gain of the reference allele.
   * @param includeGainOfReference true iff gain of reference somatic calls should be output
   * @return this, for chaining
   */
  public SomaticParamsBuilder includeGainOfReference(boolean includeGainOfReference) {
    mIncludeGainOfReference = includeGainOfReference;
    return self();
  }

  /**
   * Set site specific somatic priors.
   * @param priors the reference ranges specifying site specific somatic priors
   * @return this builder, so calls can be chained.
   */
  public SomaticParamsBuilder siteSpecificSomaticPriors(final ReferenceRanges<Double> priors) {
    mSiteSpecificSomaticPriors = priors;
    return self();
  }

  /**
   * The prior probability that a region in cancer has lost heterozygosity.
   * @param p the priors
   * @return this, for chaining
   */
  public SomaticParamsBuilder lohPrior(final double p) {
    if (p < 0 || p > 1) {
      throw new IllegalArgumentException();
    }
    mLohPrior = p;
    return self();
  }

  /**
   * The number of examples used for contamination estimation.
   * @param basis number of examples
   * @return this, for chaining
   */
  public SomaticParamsBuilder contaminationBasis(final int basis) {
    if (basis < 0) {
      throw new IllegalArgumentException();
    }
    mContaminationBasis = basis;
    return self();
  }

  /**
   * Creates a SomaticParams using the current builder configuration.
   * @return the new VariantParams
   */
  public SomaticParams create() {
    return new SomaticParams(this);
  }

  protected SomaticParamsBuilder self() {
    return this;
  }
}
