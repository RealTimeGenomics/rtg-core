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

package com.rtg.variant.bayes.snp;

import com.rtg.reference.Ploidy;
import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.complex.DescriptionComplex;
import com.rtg.variant.util.arithmetic.LogApproximatePossibility;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;

/**
 * Hypotheses for a situation where there are no hypotheses.
 */
public final class HypothesesNone<D extends Description> extends HypothesesPrior<D> {

  /** Represent a case where there are no hypotheses. */
  public static final HypothesesNone<Description> SINGLETON = new HypothesesNone<>(DescriptionNone.SINGLETON, LogApproximatePossibility.SINGLETON, 0);

  /** Represent a case where there are no hypotheses. */
  public static final HypothesesNone<DescriptionComplex> SINGLETON_COMPLEX = new HypothesesNone<>(DescriptionNone.SINGLETON_COMPLEX, LogApproximatePossibility.SINGLETON, 0);

  /**
   * @param description of the hypotheses.
   * @param arithmetic possibility arithmetic used for priors.
   * @param ref hypothesis which is the same as the reference.
   */
  public HypothesesNone(D description, PossibilityArithmetic arithmetic, int ref) {
    super(description, arithmetic, false, ref);
  }


  @Override
  public Ploidy ploidy() {
    return Ploidy.NONE;
  }
}
