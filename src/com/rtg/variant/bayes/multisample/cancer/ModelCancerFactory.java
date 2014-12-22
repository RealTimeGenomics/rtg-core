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

import com.rtg.variant.GenomePriorParams;
import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.Hypotheses;
import com.rtg.variant.bayes.ModelInterface;
import com.rtg.variant.bayes.snp.DescriptionSnp;
import com.rtg.variant.bayes.snp.HypothesesSnp;
import com.rtg.variant.bayes.snp.ModelCommonFactory;
import com.rtg.variant.bayes.snp.StatisticsSnp;
import com.rtg.variant.util.arithmetic.LogApproximatePossibility;
import com.rtg.variant.util.arithmetic.SimplePossibility;

/**
 */
public class ModelCancerFactory extends ModelCommonFactory<Description, HypothesesCancer> {

  private final double mContamination;

  /**
   * @param params information about genome used to compute priors.
   * @param contamination contamination rate.
   * @param haploid use a haploid set of hypotheses iff true.
   */
  public ModelCancerFactory(final GenomePriorParams params, final double contamination, final boolean haploid) {
    super();
    mContamination = contamination;
    for (int i = 0; i < DescriptionSnp.SINGLETON.size(); i++) {
      final HypothesesSnp hyp = new HypothesesSnp(SimplePossibility.SINGLETON, params, haploid, i);
      mHypothesesCache.add(new HypothesesCancer(hyp, LogApproximatePossibility.SINGLETON));
    }
  }

  @Override
  protected ModelInterface<Description> makeModel(final Hypotheses<Description> hyp) {
    return new ModelCancerContamination((HypothesesCancer) hyp, mContamination, new StatisticsSnp(((HypothesesCancer) hyp).subHypotheses().description()));
  }

}
