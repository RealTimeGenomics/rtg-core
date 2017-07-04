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

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.variant.bayes.Description;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;

/**
 * Provides hypotheses with priors initialized from {@code }GenomePriorParams}.
 * @param <D> class of description.
 */
@TestClass("com.rtg.variant.bayes.snp.HypothesesSnpTest")
public class HypothesesCommon<D extends Description> extends HypothesesPrior<D> {

  /**
   * @param description of the haploid
   * @param arithmetic used for reporting priors.
   * @param haploid true iff to be haploid.
   * @param ref the reference hypothesis.
   */
  public HypothesesCommon(D description, PossibilityArithmetic arithmetic, boolean haploid, int ref) {
    super(description, arithmetic, haploid, ref);
  }
}
