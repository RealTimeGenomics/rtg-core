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

package com.rtg.variant.bayes.multisample.family;

import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.Hypotheses;
import com.rtg.variant.util.arithmetic.LogApproximatePossibility;

/**
 */
public class HypothesesDisease extends Hypotheses<Description> {

  private final int mReference;

  private final double mNoDisease;

  private final double mDisease;

  /**
   * @param noDiseasePrior prior probability that there is no disease explanation.
   * @param ref the reference in a code for <code>description</code>.
   */
  HypothesesDisease(final Description description, final double noDiseasePrior, final int ref) {
    super(new DescriptionDisease(description), LogApproximatePossibility.SINGLETON, true);
    mReference = ref + 1;
    mNoDisease = arithmetic().prob2Poss(noDiseasePrior);
    mDisease = arithmetic().prob2Poss((1.0 - noDiseasePrior) / (description.size() - 1));
  }

  @Override
  public int reference() {
    return mReference;
  }

  /**
   * Gets the prior for an hypothesis.
   * @param hyp selects the hypothesis.
   * @return the prior selected by <code>hyp</code>, this will be in the format of <code>arithmetic()</code>.
   */
  public double prior(int hyp) {
    if (hyp == mReference) {
      return arithmetic().zero();
    }
    if (hyp == 0) {
      return mNoDisease;
    }
    return mDisease;
  }
}
