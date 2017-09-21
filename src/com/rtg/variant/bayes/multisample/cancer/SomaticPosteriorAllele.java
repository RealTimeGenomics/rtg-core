/*
 * Copyright (c) 2017. Real Time Genomics Limited.
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
