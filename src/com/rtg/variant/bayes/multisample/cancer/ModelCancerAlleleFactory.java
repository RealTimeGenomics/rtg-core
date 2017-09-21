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

import com.rtg.util.integrity.Exam;
import com.rtg.variant.GenomePriorParams;
import com.rtg.variant.bayes.AlleleBalanceProbability;
import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.Hypotheses;
import com.rtg.variant.bayes.HypothesesPowerSet;
import com.rtg.variant.bayes.ModelInterface;
import com.rtg.variant.bayes.snp.DescriptionSnp;
import com.rtg.variant.bayes.snp.HypothesesSnp;
import com.rtg.variant.bayes.snp.ModelCommonFactory;
import com.rtg.variant.bayes.snp.StatisticsSnp;
import com.rtg.variant.util.arithmetic.LogApproximatePossibility;

/**
 * A factory capable of generating model instances for cancer calling.
 */
public class ModelCancerAlleleFactory extends ModelCommonFactory<Description, HypothesesPowerSet<Description>> {

  /**
   * @param params information about genome used to compute priors.
   * @param haploid use a haploid set of hypotheses iff true.
   * @param alleleBalance allele balance probability implementation
   */
  public ModelCancerAlleleFactory(final GenomePriorParams params, final boolean haploid, final AlleleBalanceProbability alleleBalance) {
    super(alleleBalance);
    final HypothesesSnp unknownHypothesesSnp = new HypothesesSnp(LogApproximatePossibility.SINGLETON, params, haploid, Hypotheses.NO_HYPOTHESIS);
    mHypothesisUnknown = new HypothesesPowerSet<>(unknownHypothesesSnp.description(), LogApproximatePossibility.SINGLETON, unknownHypothesesSnp.reference());
    for (int i = 0; i < DescriptionSnp.SINGLETON.size(); ++i) {
      //final HypothesesSnp hyp = new HypothesesSnp(SimplePossibility.SINGLETON, params, haploid, i);
      mHypothesesCache.add(new HypothesesPowerSet<>(DescriptionSnp.SINGLETON, LogApproximatePossibility.SINGLETON, i));
    }
  }

  @Override
  protected ModelInterface<Description> makeModel(final Hypotheses<Description> hyp) {
    assert hyp instanceof HypothesesPowerSet;
    final HypothesesPowerSet<Description> hypothesesAllele = (HypothesesPowerSet<Description>) hyp;
    return new ModelCancerAllele<>(hypothesesAllele, new StatisticsSnp(hypothesesAllele.description()));
  }

  @Override
  public boolean globalIntegrity() {
    for (int i = 0; i < mHypothesesCache.size(); ++i) {
      Exam.assertEquals((1 << i) - 1, mHypothesesCache.get(i).reference());
    }
    return integrity();
  }
}
