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
