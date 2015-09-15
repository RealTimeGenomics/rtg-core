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
class PosteriorPure extends AbstractPosterior {

  /**
   * @param qa The Q matrix of cancer mutation probabilities.
   * @param normal array of normal hypotheses
   * @param cancer array of cancer hypotheses
   * @param hypotheses the hypotheses containing priors
   */
  PosteriorPure(final double[][] qa, final ModelInterface<?> normal, final ModelInterface<?> cancer, HypothesesPrior<?> hypotheses) {
    super(hypotheses);
    //System.err.println("normal " + normal);
    //System.err.println("cancer " + cancer);
    for (int i = 0; i < mLength; i++) {
      final double pi = hypotheses.arithmetic().poss2Ln(hypotheses.p(i)) + normal.posteriorLn0(i);
      for (int j = 0; j < mLength; j++) {
        final double pj = cancer.posteriorLn0(j);
        final double q = MathUtils.log(qa[i][j]);
        final double t = q + pi + pj;
        mPosterior[i][j] = t;
        //System.err.println("PosteriorPure i=" + i + " j=" + j + " hypNormal=" + hypotheses.name(i) + " hypCancer=" + hypotheses.name(j) + " q=" + Utils.realFormat(q, 3) + " pi=" + Utils.realFormat(pi, 3) + " pj=" + Utils.realFormat(pj, 3) + " t=" + Utils.realFormat(t, 3));
      }
    }
    final double r = Math.log(0.0001); // XXX roughly should be prob of machine error -- via constructor?
    contraryEvidenceAdjustment(normal.statistics(), cancer.statistics(), r, r);
    postConstruction();
  }

}
