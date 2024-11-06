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

package com.rtg.variant.bayes.multisample.cancer;

import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.Hypotheses;
import com.rtg.variant.bayes.snp.DescriptionCommon;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;

/**
 * Asymmetric pairs of hypotheses used in contaminated cancer model.
 * @param <S> type of the underlying hypotheses
 */
public class HypothesesCancer<S extends Hypotheses<? extends Description>> extends Hypotheses<Description> {

  private final S mSubHypotheses;

  /**
   * @param hypotheses underlying hypotheses pairs of which form the cancer hypotheses.
   * @param arithmetic used in calculations and for priors.
   */
  protected HypothesesCancer(S hypotheses, PossibilityArithmetic arithmetic) {
    super(makeCancerDescription(hypotheses), arithmetic, new CodeCross(hypotheses.size()), false);
    mSubHypotheses = hypotheses;
  }

  private static Description makeCancerDescription(Hypotheses<? extends Description> hypotheses) {
    final String[] names = new String[hypotheses.size()];
    for (int i = 0; i < hypotheses.size(); ++i) {
      names[i] = hypotheses.name(i);
    }
    return new DescriptionCommon(names);
  }

  @Override
  public int reference() {
    return subHypotheses().reference();
  }

  /**
   * @return the underlying hypotheses.
   */
  public S subHypotheses() {
    return mSubHypotheses;
  }

}
