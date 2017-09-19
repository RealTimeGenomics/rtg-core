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

import com.rtg.util.MathUtils;
import com.rtg.variant.bayes.ModelInterface;
import com.rtg.variant.bayes.snp.HypothesesPrior;

/**
 * Calculate the joint posterior and marginal distributions for cancer.
 */
class SomaticPosteriorPure extends AbstractSomaticPosterior {

  /**
   * @param q The Q matrix of cancer mutation probabilities.
   * @param normal array of normal hypotheses
   * @param cancer array of cancer hypotheses
   * @param hypotheses the hypotheses containing priors
   * @param phi probability of seeing contrary evidence in the original
   * @param psi probability of seeing contrary evidence in the derived
   */
  SomaticPosteriorPure(final double[][] q, final ModelInterface<?> normal, final ModelInterface<?> cancer, HypothesesPrior<?> hypotheses, double phi, double psi) {
    super(hypotheses, phi, psi);
    //System.err.println("normal " + normal);
    //System.err.println("cancer " + cancer);
    for (int normalHyp = 0; normalHyp < q.length; ++normalHyp) {
      final double pNormal = hypotheses.arithmetic().poss2Ln(hypotheses.p(normalHyp)) + normal.posteriorLn0(normalHyp);
      for (int cancerHyp = 0; cancerHyp < q[normalHyp].length; ++cancerHyp) {
        final double pCancer = cancer.posteriorLn0(cancerHyp);
        final double qv = MathUtils.log(q[normalHyp][cancerHyp]);
        final double t = qv + pNormal + pCancer;
        mPosterior[normalHyp][cancerHyp] = t;
      }
    }
    contraryEvidenceAdjustment(normal.statistics(), cancer.statistics());
    postConstruction();
  }

}
