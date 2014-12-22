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
import com.rtg.variant.bayes.Description;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;

/**
 * Encapsulates the space of SNP hypotheses
 */
public final class HypothesesSnp extends HypothesesCommon<Description> {

  /**
   * @param arithmetic possibility arithmetic used for priors.
   * @param params parameters used for setting priors.
   * @param haploid true iff the hypotheses are to be haploid (rather than diploid).
   * @param ref hypothesis which is the same as the reference.
   */
  public HypothesesSnp(PossibilityArithmetic arithmetic, final GenomePriorParams params, final boolean haploid, int ref) {
    super(DescriptionSnp.SINGLETON, arithmetic, haploid, ref);
    if (ref == -1) {
      initPriors(params);
    } else {
      initPriors(params, ref);
    }
  }

  private void initPriors(final GenomePriorParams params, int ref) {
    double total = 0.0;
    for (int i = 0; i < size(); i++) {
      final String name = name(i);
      final double prob = params.getPriorDistr(name)[ref];
      total += prob;
    }
    for (int i = 0; i < size(); i++) {
      final String name = name(i);
      final double prob = params.getPriorDistr(name)[ref] / total;
      final double prior = arithmetic().prob2Poss(prob);
      setPrior(i, prior);
    }
  }

  private void initPriors(final GenomePriorParams params) {
    double total = 0.0;
    for (int i = 0; i < size(); i++) {
      final String name = name(i);
      for (int ref = 0; ref < DescriptionSnp.SINGLETON.size(); ref++) {
        final double prob = params.getPriorDistr(name)[ref];
        total += prob;
      }
    }
    for (int i = 0; i < size(); i++) {
      final String name = name(i);
      double prob = 0.0;
      for (int ref = 0; ref < DescriptionSnp.SINGLETON.size(); ref++) {
        prob += params.getPriorDistr(name)[ref];
      }
      prob = prob / total;
      final double prior = arithmetic().prob2Poss(prob);
      setPrior(i, prior);
    }
  }

}
