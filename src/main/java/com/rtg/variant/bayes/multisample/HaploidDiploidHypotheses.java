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

package com.rtg.variant.bayes.multisample;

import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.ModelInterface;
import com.rtg.variant.bayes.multisample.population.DescriptionCounts;
import com.rtg.variant.bayes.snp.HypothesesPrior;

/**
 * A wrapper class that holds both a haploid and diploid set of hypotheses.
 * @param <T> type of hypotheses
 */
public class HaploidDiploidHypotheses<T extends HypothesesPrior<? extends Description>> {

  public DescriptionCounts getDescriptionCounts() {
    return mDescriptionCounts;
  }

  private final DescriptionCounts mDescriptionCounts;

  private final T mHaploid;
  private final T mDiploid;
  private final boolean mIsDefault;
  private final T mNone;

  /**
   * Constructor for hypotheses that originate from default genome priors.
   * Will compute new priors.
   * @param none set of empty hypotheses.
   * @param haploid set of haploid hypotheses
   * @param diploid set of diploid hypotheses
  */
  public HaploidDiploidHypotheses(T none, T haploid, T diploid) {
    this(none, haploid, diploid, true, null);
  }

  /**
   * @param none set of empty hypotheses.
   * @param haploid set of haploid hypotheses
   * @param diploid set of diploid hypotheses
   * @param isDefault true if these hypotheses are come from default genome priors without optimization
   * @param descriptionCounts description counts object
   */
  public HaploidDiploidHypotheses(T none, T haploid, T diploid, boolean isDefault, DescriptionCounts descriptionCounts) {
    if (haploid != null && !haploid.haploid() || diploid != null && diploid.haploid()) {
      throw new IllegalArgumentException();
    }
    mNone = none;
    mHaploid = haploid;
    mDiploid = diploid;
    mIsDefault = isDefault;
    mDescriptionCounts = descriptionCounts;
    assert mHaploid == null || mDescriptionCounts == null || mDescriptionCounts.getSize() >= mHaploid.size();
  }

  /**
   * Get whether the priors were created from default genome priors.
   * @return Returns true if the priors were created from default genome priors.
   */
  public boolean isDefault() {
    return mIsDefault;
  }

  /**
   * Get the set of hypotheses appropriate for the model.
   * @param model the model
   * @return the hypotheses selected based on hypotheses contained within the model
   */
  public T get(ModelInterface<?> model) {
    switch (model.hypotheses().ploidy()) {
      case NONE:
        return mNone;
      case HAPLOID:
        return mHaploid;
      case DIPLOID:
        return mDiploid;
      case POLYPLOID:
      default:
        throw new RuntimeException();
    }
    //return model.hypotheses() == HypothesesNone.SINGLETON ? HypothesesNone.SINGLETON : model.haploid() ? mHaploid : mDiploid;
  }

  /**
   *
   * @return the haploid set of hypotheses
   */
  public T haploid() {
    return mHaploid;
  }

  /**
   * @return the diploid set of hypotheses
   */
  public T diploid() {
    return mDiploid;
  }

  /**
   * @return the set of empty of hypotheses
   */
  public T none() {
    return mNone;
  }
}
