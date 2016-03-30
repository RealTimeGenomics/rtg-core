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
import com.rtg.variant.bayes.AlleleBalanceProbability;
import com.rtg.variant.bayes.Code;
import com.rtg.variant.bayes.ModelInterface;
import com.rtg.variant.bayes.snp.HypothesesPrior;

/**
 * Calculate the joint posterior and marginal distributions for cancer.
 */
class SomaticPosteriorContaminated extends AbstractSomaticPosterior {

  /**
   * @param qa The Q matrix of cancer mutation probabilities.
   * @param normal model.
   * @param cancer model, has cross product hypotheses to allow for contamination.
   * @param hypotheses the hypotheses containing priors
   * @param phi probability of seeing contrary evidence in the original
   * @param psi probability of seeing contrary evidence in the derived
   */
  SomaticPosteriorContaminated(final double[][] qa, final ModelInterface<?> normal, final ModelInterface<?> cancer, HypothesesPrior<?> hypotheses, double phi, double psi) {
    super(normal.hypotheses(), phi, psi);
    //System.err.println("normal " + normal);
    //System.err.println("cancer " + cancer);
    assert cancer instanceof ModelCancerContamination;
    assert !cancer.haploid(); // The cancer model is a cross-product of normal x cancer hypotheses
    final Code code = cancer.hypotheses().code();
    final AlleleBalanceProbability alleleBalance = normal.alleleBalanceProbability();
    for (int i = 0; i < mPosterior.length; i++) {
      final double pi = hypotheses.arithmetic().poss2Ln(hypotheses.p(i)) + normal.posteriorLn0(i);
      for (int j = 0; j < mPosterior.length; j++) {
        final int k = code.code(i, j); // code point for normal(i) x cancer(j)
        assert k >= 0 && k < code.size() : k + " " + code.size();
        final double pj = cancer.posteriorLn0(k);
        final double q = MathUtils.log(qa[i][j]) + mArithmetic.ln2Poss(alleleBalance.alleleBalanceLn(i, normal.hypotheses(), cancer.statistics()));
        final double t = q + pi + pj;
        mPosterior[i][j] = t;
      }
    }
    // After this point the special code for contamination no longer plays any part, mPosterior
    // is normal x cancer for both contaminated and uncontaminated callers.
    contraryEvidenceAdjustment(normal.statistics(), cancer.statistics());
    postConstruction();
  }

}
