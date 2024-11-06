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

import java.util.List;

import com.rtg.mode.DnaUtils;
import com.rtg.reference.Ploidy;
import com.rtg.variant.Variant;
import com.rtg.variant.VariantLocus;
import com.rtg.variant.VariantOutputLevel;
import com.rtg.variant.VariantOutputOptions;
import com.rtg.variant.VariantParams;
import com.rtg.variant.VariantSample;
import com.rtg.variant.bayes.AlleleStatistics;
import com.rtg.variant.bayes.Code;
import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.GenotypeMeasure;
import com.rtg.variant.bayes.Hypotheses;
import com.rtg.variant.bayes.ModelInterface;
import com.rtg.variant.bayes.multisample.HaploidDiploidHypotheses;
import com.rtg.variant.bayes.snp.HypothesesPrior;
import com.rtg.mode.DNARange;
import com.rtg.variant.util.VariantUtils;

/**
 * An allele based somatic caller.
 */
public class AlleleSomaticCaller extends AbstractSomaticCaller {

  /**
   * An allele based somatic caller.
   * @param qHaploidFactory Q matrix factory for haploid normal
   * @param qDiploidFactory Q matrix factory for diploid normal
   * @param params variant params
   * @param phi probability of seeing contrary evidence in the normal
   * @param psi probability of seeing contrary evidence in the cancer
   */
  public AlleleSomaticCaller(final AlleleSomaticPriorsFactory<?> qHaploidFactory, final AlleleSomaticPriorsFactory<?> qDiploidFactory, final VariantParams params, final double phi, final double psi) {
    super(qHaploidFactory, qDiploidFactory, params, phi, psi);
  }

  @Override
  protected AbstractSomaticPosterior makePosterior(final ModelInterface<?> normal, final ModelInterface<?> cancer, final HypothesesPrior<?> normalPrior, final double mu) {
    return new SomaticPosteriorAllele((normal.haploid() ? mQHaploidFactory : mQDiploidFactory).somaticQ(mu), normal, cancer, normalPrior, mPhi, mPsi);
  }

  /* Convert from normal hypothesis to corresponding cancer hypothesis. */
  private int toCancerHyp(final Code normalCode, final int normalHyp) {
    final int a = normalCode.a(normalHyp);
    final int b = normalCode.bc(normalHyp);
    return ((1 << a) | (1 << b)) - 1; // Depends on details of HypothesesPowerSet
  }

  private Ploidy cancerPloidy(final Ploidy normalPloidy, final int cancerHyp) {
    final int alleleCount = Integer.bitCount(cancerHyp + 1); // Depends on HypothesesPowerSet coding
    if (alleleCount > 2) {
      return Ploidy.POLYPLOID;
    }
    if (alleleCount == 2) {
      return Ploidy.DIPLOID;
    }
    // We can have a cancer call with just one allele.  Although such a situation could correspond
    // to a copy number variation or loss of heterozygosity, it is probably best to make these
    // calls come out with the same ploidy as the normal sample at this site.
    return normalPloidy;
  }

  private VariantSample setCancerCallValues(final GenotypeMeasure posterior, final int cancerHyp, final ModelInterface<?> cancerModel, final VariantOutputOptions params, final Ploidy normalPloidy, final VariantSample.DeNovoStatus dns, final Double dnp) {
    final Hypotheses<?> cancerHypotheses = cancerModel.hypotheses();
    final Ploidy ploidy = cancerPloidy(normalPloidy, cancerHyp);
    String cancerHypName = cancerHypotheses.name(cancerHyp);
    // When making a homozygous diploid call, we need to double up the allele
    if (ploidy == Ploidy.DIPLOID && cancerHypName.indexOf(VariantUtils.COLON) < 0) {
      cancerHypName = cancerHypName + VariantUtils.COLON + cancerHypName;
    }
    final VariantSample sample = new VariantSample(ploidy, cancerHypName, cancerHypotheses.reference() == cancerHyp, posterior, dns, dnp);
    cancerModel.statistics().addCountsToSample(sample, cancerModel, params);
    return sample;
  }

  @Override
  public <D extends Description, T extends HypothesesPrior<D>> Variant makeCall(String templateName, int position, int endPosition, byte[] ref, List<ModelInterface<?>> models, HaploidDiploidHypotheses<T> normalHypotheses) {
    assert DNARange.valid(ref[position], DNARange.N, DNARange.T);
    assert models.size() == 2;

    final ModelInterface<?> modelNormal = models.get(NORMAL);
    final ModelInterface<?> modelCancer = models.get(CANCER);
    assert modelCancer instanceof ModelCancerAllele;
    final HypothesesPrior<?> hypotheses = normalHypotheses.get(modelNormal);
    final Code code = hypotheses.code();

    final AbstractSomaticPosterior posterior = makePosterior(modelNormal, modelCancer, hypotheses, getSomaticPrior(templateName, position, endPosition));
    final boolean isSomatic;
    final int bestNormalHyp = posterior.bestNormal();
    final int bestCancerHyp = posterior.bestCancer();
    final boolean sameCall = posterior.isSameCall() || (toCancerHyp(code, bestNormalHyp) == bestCancerHyp);
    final Ploidy normalPloidy = hypotheses.haploid() ? Ploidy.HAPLOID : Ploidy.DIPLOID;
    final String refAllele = DnaUtils.bytesToSequenceIncCG(ref, position, endPosition - position);
    final double ratio = posterior.posteriorScore();

    final AlleleStatistics<?> ac = modelCancer.statistics().counts();
    final Description description = ac.getDescription();
    final int va = mVariantAlleleTrigger.getVariantAllele(ac, description, refAllele);

    final boolean triggersVa = va != -1;
    if (sameCall && hypotheses.reference() == bestNormalHyp && ratio >= mIdentityInterestingThreshold  // Call was same for both samples and equal to the reference
      && !(triggersVa || (mParams.callLevel() == VariantOutputLevel.ALL))) {
      // We don't need to output any record here.
      return null;
    }

    boolean interesting = true;
    if (sameCall) {
      if (hypotheses.reference() == bestNormalHyp && ratio >= mIdentityInterestingThreshold) {
        interesting = false;
      }
      isSomatic = false;
    } else {
      isSomatic = true;
    }

    final VariantSample normalSample = setCallValues(posterior.normalMeasure(), bestNormalHyp, hypotheses, modelNormal, mParams, normalPloidy, VariantSample.DeNovoStatus.UNSPECIFIED, null);
    final VariantSample cancerSample;
    if (isSomatic) {
      cancerSample = setCancerCallValues(posterior.cancerMeasure(), bestCancerHyp, modelCancer, mParams, normalPloidy, VariantSample.DeNovoStatus.IS_DE_NOVO, ratio);
    } else {
      cancerSample = setCancerCallValues(posterior.cancerMeasure(), bestCancerHyp, modelCancer, mParams, normalPloidy, VariantSample.DeNovoStatus.NOT_DE_NOVO, null);
    }
    if (triggersVa) {
      cancerSample.setVariantAllele(description.name(va));
    }
    final VariantLocus locus = new VariantLocus(templateName, position, endPosition, refAllele, VariantUtils.getPreviousRefNt(ref, position));
    final Variant v = new Variant(locus, normalSample, cancerSample);
    if (modelNormal.statistics().overCoverage(mParams, templateName) || modelCancer.statistics().overCoverage(mParams, templateName)) {
      v.addFilter(Variant.VariantFilter.COVERAGE);
    } else if (modelNormal.statistics().ambiguous(mParams) || modelCancer.statistics().ambiguous(mParams)) {
      v.addFilter(Variant.VariantFilter.AMBIGUITY);
    }
    if (interesting) {
      v.setInteresting();
    }
    if (isSomatic || mParams.somaticParams().includeGermlineVariants()) {
      v.setNormalCancerScore(posterior.ncScore());
    }
    //updateParameterEstimationCounts(normalPloidy, code, bestNormal, bestCancer, loh); // todo do we really care -- its likely bogus anyway
    return v;
  }

}
