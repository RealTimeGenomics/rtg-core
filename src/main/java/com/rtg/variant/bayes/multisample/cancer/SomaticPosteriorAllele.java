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

import com.rtg.util.MathUtils;
import com.rtg.variant.bayes.ModelInterface;
import com.rtg.variant.bayes.snp.HypothesesPrior;

/**
 * Posterior calculations for allele based cancer calling.
 */
public class SomaticPosteriorAllele extends AbstractSomaticPosterior {

  /**
   * Posterior for the allele based cancer caller.
   * @param q The Q matrix of cancer mutation probabilities.
   * @param normal normal sample model
   * @param cancer cancer sample model
   * @param normalPrior priors for normal sample
   * @param phi probability of seeing contrary evidence in the normal
   * @param psi probability of seeing contrary evidence in the cancer
   */
  public SomaticPosteriorAllele(final double[][] q, final ModelInterface<?> normal, final ModelInterface<?> cancer, final HypothesesPrior<?> normalPrior, final double phi, final double psi) {
    super(normal.hypotheses(), cancer.hypotheses(), phi, psi);
    for (int normalHyp = 0; normalHyp < normal.size(); ++normalHyp) {
      final double pNormal = normal.arithmetic().poss2Ln(normalPrior.p(normalHyp)) + normal.posteriorLn0(normalHyp);
      for (int cancerHyp = 0; cancerHyp < cancer.size(); ++cancerHyp) {
        final double pCancer = cancer.posteriorLn0(cancerHyp);
        final double qv = MathUtils.log(q[normalHyp][cancerHyp]);
        mPosterior[normalHyp][cancerHyp] = qv + pNormal + pCancer;
      }
    }
    contraryEvidenceAdjustment(normal.statistics(), cancer.statistics());
    postConstruction();
  }

  @Override
  protected boolean isSame(final int normalHyp, final int cancerHyp) {
    return normalHypToAlleleBits(normalHyp) == cancerHypToAlleleBits(cancerHyp);
  }

  @Override
  protected int cancerHypToAlleleBits(final int cancerHyp) {
    return cancerHyp + 1;
  }
}
