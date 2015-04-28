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
import com.rtg.variant.bayes.Code;
import com.rtg.variant.bayes.ModelInterface;
import com.rtg.variant.bayes.snp.HypothesesPrior;

/**
 * Calculate the joint posterior and marginal distributions for cancer.
 */
class PosteriorContaminated extends AbstractPosterior {

  /**
   * @param qa The Q matrix of cancer mutation probabilities.
   * @param normal model.
   * @param cancer model, has cross product hypotheses to allow for contamination.
   * @param hypotheses the hypotheses containing priors
   */
  PosteriorContaminated(final double[][] qa, final ModelInterface<?> normal, final ModelInterface<?> cancer, HypothesesPrior<?> hypotheses) {
    super(normal.hypotheses());
    //System.err.println("normal " + normal);
    //System.err.println("cancer " + cancer);
    assert !cancer.haploid();
    final Code code = cancer.hypotheses().code();
    for (int i = 0; i < mPosterior.length; i++) {
      final double pi = hypotheses.arithmetic().poss2Ln(hypotheses.p(i)) + normal.posteriorLn0(i);
      for (int j = 0; j < mPosterior.length; j++) {
        final int k = code.code(i, j);
        assert k >= 0 && k < code.size() : k + " " + code.size();
        final double pj = cancer.posteriorLn0(k);
        final double q = MathUtils.log(qa[i][j]);
        final double t = q + pi + pj;
        //System.err.println("PosteriorContaminated i=" + i + " j=" + j + " hypNormal=" + hypotheses.name(i) + " hypCancer=" + hypotheses.name(j) + " q=" + Utils.realFormat(q, 3) + " pi=" + Utils.realFormat(pi, 3) + " pj=" + Utils.realFormat(pj, 3) + " t=" + Utils.realFormat(t, 3));
        mPosterior[i][j] = t;
      }
    }
    postConstruction();
  }

}
