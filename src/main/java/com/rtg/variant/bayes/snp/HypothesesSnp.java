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

package com.rtg.variant.bayes.snp;

import com.rtg.variant.GenomePriorParams;
import com.rtg.variant.bayes.Description;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;

/**
 * Encapsulates the space of SNP hypotheses.
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
    if (ref == NO_HYPOTHESIS) {
      initPriors(params);
    } else {
      initPriors(params, ref);
    }
  }

  private void initPriors(final GenomePriorParams params, int ref) {
    double total = 0.0;
    for (int i = 0; i < size(); ++i) {
      final String name = name(i);
      final double prob = params.getPriorDistr(name)[ref];
      total += prob;
    }
    for (int i = 0; i < size(); ++i) {
      final String name = name(i);
      final double prob = params.getPriorDistr(name)[ref] / total;
      final double prior = arithmetic().prob2Poss(prob);
      setPrior(i, prior);
    }
  }

  private void initPriors(final GenomePriorParams params) {
    double total = 0.0;
    for (int i = 0; i < size(); ++i) {
      final String name = name(i);
      for (int ref = 0; ref < DescriptionSnp.SINGLETON.size(); ++ref) {
        final double prob = params.getPriorDistr(name)[ref];
        total += prob;
      }
    }
    for (int i = 0; i < size(); ++i) {
      final String name = name(i);
      double prob = 0.0;
      for (int ref = 0; ref < DescriptionSnp.SINGLETON.size(); ++ref) {
        prob += params.getPriorDistr(name)[ref];
      }
      prob = prob / total;
      final double prior = arithmetic().prob2Poss(prob);
      setPrior(i, prior);
    }
  }

}
