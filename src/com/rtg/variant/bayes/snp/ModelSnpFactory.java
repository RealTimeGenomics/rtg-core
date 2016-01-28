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

import com.rtg.variant.GenomePriorParams;
import com.rtg.variant.bayes.AlleleBalanceProbability;
import com.rtg.variant.bayes.Description;
import com.rtg.variant.util.arithmetic.LogApproximatePossibility;

/**
 */
public class ModelSnpFactory extends ModelCommonFactory<Description, HypothesesSnp> {

  /**
   * @param params information about genome used to compute priors.
   * @param haploid use a haploid set of hypotheses iff true.
   * @param alleleBalance allele balance probability implementation
   */
  public ModelSnpFactory(final GenomePriorParams params, final boolean haploid, final AlleleBalanceProbability alleleBalance) {
    super(alleleBalance);
    mHypothesisUnknown = new HypothesesSnp(LogApproximatePossibility.SINGLETON, params, haploid, -1);
    for (int i = 0; i < DescriptionSnp.SINGLETON.size(); i++) {
      mHypothesesCache.add(new HypothesesSnp(LogApproximatePossibility.SINGLETON, params, haploid, i));
    }
  }
}
