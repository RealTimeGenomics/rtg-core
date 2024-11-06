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

package com.rtg.variant.bayes.multisample.family;

import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.Hypotheses;
import com.rtg.variant.util.arithmetic.LogApproximatePossibility;

/**
 * Set of hypotheses with a disease hypothesis.
 */
public class HypothesesDisease extends Hypotheses<Description> {

  private final int mReference;

  private final double mNoDisease;

  private final double mDisease;

  /**
   * @param description the base descriptions
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
