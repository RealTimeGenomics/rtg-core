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
