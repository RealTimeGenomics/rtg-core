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
    for (int i = 0; i < hypotheses.size(); i++) {
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
