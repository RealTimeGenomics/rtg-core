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

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.mode.DnaUtils;
import com.rtg.reference.Ploidy;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;
import com.rtg.util.intervals.RangeList;
import com.rtg.util.intervals.ReferenceRanges;
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
import com.rtg.variant.bayes.multisample.MultisampleJointCaller;
import com.rtg.variant.bayes.multisample.VariantAlleleTrigger;
import com.rtg.variant.bayes.snp.HypothesesPrior;
import com.rtg.mode.DNARange;
import com.rtg.variant.util.VariantUtils;

/**
 * The bulk of the implementation for the somatic caller.  Takes the models from the normal and cancer
 * samples and examines the joint (possibly contaminated) distribution to determine the best call for
 * the pair of samples.
 */
@TestClass(value = {"com.rtg.variant.bayes.multisample.cancer.ContaminatedSomaticCallerTest", "com.rtg.variant.bayes.multisample.cancer.PureSomaticCallerTest"})
public abstract class AbstractSomaticCaller extends IntegralAbstract implements MultisampleJointCaller {

  static final int NORMAL = 0;
  static final int CANCER = 1;

  protected final SomaticPriorsFactory mQHaploidFactory;
  protected final SomaticPriorsFactory mQDiploidFactory;
  protected final VariantParams mParams;
  private final ReferenceRanges<Double> mSiteSpecificSomaticPriors;
  protected final double mIdentityInterestingThreshold;
  protected final double mPhi;
  protected final double mPsi;
  protected final VariantAlleleTrigger mVariantAlleleTrigger;

  /**
   * @param qHaploidFactory haploid Q matrix factory
   * @param qDiploidFactory diploid Q matrix factory
   * @param params variant params
   * @param phi probability of seeing contrary evidence in the original
   * @param psi probability of seeing contrary evidence in the derived
   */
  public AbstractSomaticCaller(final SomaticPriorsFactory qHaploidFactory, final SomaticPriorsFactory qDiploidFactory, final VariantParams params, final double phi, final double psi) {
    mQHaploidFactory = qHaploidFactory;
    mQDiploidFactory = qDiploidFactory;
    mParams = params;
    mSiteSpecificSomaticPriors = mParams.somaticParams().siteSpecificSomaticPriors();
    mIdentityInterestingThreshold = mParams.interestingThreshold();
    mPhi = phi;
    mPsi = psi;
    mVariantAlleleTrigger = new VariantAlleleTrigger(params.minVariantAllelicDepth(), params.minVariantAllelicFraction());
  }

  /**
   * Construct an appropriate posterior. Differs in the contaminated and non-contaminated case.
   * @param normal bayesian for the normal genome
   * @param cancer bayesian for the cancer genome
   * @param hypotheses the hypotheses containing priors
   * @param mu somatic mutation rate
   * @return the posterior
   */
  protected abstract AbstractSomaticPosterior makePosterior(final ModelInterface<?> normal, final ModelInterface<?> cancer, final HypothesesPrior<?> hypotheses, final double mu);

  protected VariantSample setCallValues(GenotypeMeasure posterior, int cat, Hypotheses<?> hypotheses, ModelInterface<?> model, VariantOutputOptions params, Ploidy ploidy, VariantSample.DeNovoStatus dns, Double dnp) {
    final VariantSample sample = new VariantSample(ploidy, hypotheses.name(cat), hypotheses.reference() == cat, posterior, dns, dnp);
    model.statistics().addCountsToSample(sample, model, params);
    return sample;
  }

  private double getSomaticPrior(final String seqName, final int pos) {
    if (mSiteSpecificSomaticPriors != null) {
      final RangeList<Double> rangeList = mSiteSpecificSomaticPriors.get(seqName);
      if (rangeList != null) {
        final List<Double> v = rangeList.find(pos);
        // Take the maximum of the supplied priors
        if (v != null) {
          double p = 0;
          for (final double pv : v) {
            if (pv > p) {
              p = pv;
            }
          }
          return p;
        }
      }
    }
    return mParams.somaticParams().somaticRate();
  }

  protected double getSomaticPrior(final String seqName, final int start, final int end) {
    if (end <= start + 1) {
      return getSomaticPrior(seqName, start); // handles zero and one length regions
    }
    // For an extended region choose a prior that is the mean of the point priors in the region
    double s = 0;
    for (int k = start; k < end; ++k) {
      s += getSomaticPrior(seqName, k);
    }
    return s / (end - start);
  }

  private double loh(final Hypotheses<?> hypotheses, final int normal, final int cancer) {
    final Code code = hypotheses.code();
    if (!code.homozygous(normal) && code.homozygous(cancer)) {
      return 1;
    }
    if (cancer == normal) {
      return 0;
    }
    return -1;
  }

  // Counts for making LOH and somatic rate estimations, see multiscoring.tex theory doc
  private double mEll = 0;
  private double mOneMinusEll = 0;
  private double mMu = 0;
  private double mOneMinusMu = 0;

  private void updateParameterEstimationCounts(final Ploidy normalPloidy, final Code code, final int normal, final int cancer, final double loh) {
    // See multiscoring.tex for an explanation.
    // A more thorough treatment should probably increment mEll and mOneMinusEll in proportion to loh
    // but we often have loh of 0 when we haven't actually made a call, so this could be tricky.
    // The theory really only covers the SNP cases, but here we just "make it work" for complex
    // cases as well.
    if (loh <= 0) {
      ++mOneMinusEll;
      if (normalPloidy == Ploidy.HAPLOID) {
        if (normal == cancer) {
          ++mOneMinusMu;
        } else {
          ++mMu;
        }
      } else {
        if (normal == cancer) {
          mOneMinusMu += 2; // efficiency, this case would also be handled below
        } else {
          final int na = code.a(normal);
          final int nb = code.bc(normal);
          final int ca = code.a(cancer);
          final int cb = code.bc(cancer);
          final int m = ((na == ca || na == cb) ? 1 : 0) + ((nb == ca || nb == cb) ? 1 : 0);
          mOneMinusMu += m;
          mMu += 2 - m;
        }
      }
    } else {
      assert code.homozygous(cancer);
      ++mEll;
      if (code.a(normal) == cancer || code.bc(normal) == cancer) {
        ++mOneMinusMu;
      } else {
        ++mMu;
      }
    }
  }

  @Override
  public <D extends Description, T extends HypothesesPrior<D>> Variant makeCall(String templateName, int position, int endPosition, byte[] ref, List<ModelInterface<?>> models, HaploidDiploidHypotheses<T> normalHypotheses) {
    assert DNARange.valid(ref[position], DNARange.N, DNARange.T);
    assert models.size() == 2;

    final ModelInterface<?> modelNormal = models.get(NORMAL);
    final ModelInterface<?> modelCancer = models.get(CANCER);
    final HypothesesPrior<?> hypotheses = normalHypotheses.get(modelNormal);
    final Code code = hypotheses.code();

    final AbstractSomaticPosterior posterior = makePosterior(modelNormal, modelCancer, hypotheses, getSomaticPrior(templateName, position, endPosition));
    final boolean isSomatic;
    final int bestNormal = posterior.bestNormal();
    final int bestCancer = posterior.bestCancer();
    final boolean sameCall = posterior.isSameCall() || (bestNormal == bestCancer);
    // Simple LOH test based on ploidy of results alone, could be done with Bayesian calculation later
    final double loh = loh(hypotheses, bestNormal, bestCancer);
    final Ploidy normalPloidy = hypotheses.haploid() ? Ploidy.HAPLOID : Ploidy.DIPLOID;
    final String refAllele = DnaUtils.bytesToSequenceIncCG(ref, position, endPosition - position);
    final double ratio = posterior.posteriorScore();


    final AlleleStatistics<?> ac = modelCancer.statistics().counts();
    final Description description = ac.getDescription();
    final int va = mVariantAlleleTrigger.getVariantAllele(ac, description, refAllele);

    final boolean triggersVa = va != -1;
    if (sameCall && hypotheses.reference() == bestNormal && ratio >= mIdentityInterestingThreshold  // Call was same for both samples and equal to the reference
      && !(triggersVa || (mParams.callLevel() == VariantOutputLevel.ALL))) {
      // We don't need to output any record here.
      return null;
    }

    boolean interesting = true;
    if (sameCall) {
      if (hypotheses.reference() == bestNormal && ratio >= mIdentityInterestingThreshold) {
        interesting = false;
      }
      isSomatic = false;
    } else {
      isSomatic = true;
    }

    final VariantSample normalSample = setCallValues(posterior.normalMeasure(), bestNormal, hypotheses, modelNormal, mParams, normalPloidy, VariantSample.DeNovoStatus.UNSPECIFIED, null);
    final VariantSample cancerSample;
    if (isSomatic) {
      cancerSample = setCallValues(posterior.cancerMeasure(), bestCancer, hypotheses, modelCancer, mParams, normalPloidy, VariantSample.DeNovoStatus.IS_DE_NOVO, ratio);
    } else {
      cancerSample = setCallValues(posterior.cancerMeasure(), bestCancer, hypotheses, modelCancer, mParams, normalPloidy, VariantSample.DeNovoStatus.NOT_DE_NOVO, null);
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
    updateParameterEstimationCounts(normalPloidy, code, bestNormal, bestCancer, loh);
    return v;
  }

  @Override
  public void toString(StringBuilder sb) {
    // Note this dump of Q does not deal with any somatic site-specific priors that might be active
    final double[][] qMatrix = (mQHaploidFactory != null ? mQHaploidFactory : mQDiploidFactory).somaticQ(mParams.somaticParams().somaticRate());
    sb.append("length=").append(qMatrix.length).append(LS);
    for (final double[] q : qMatrix) {
      for (final double v : q) {
        sb.append(com.rtg.util.Utils.realFormat(v, 3)).append(" ");
      }
      sb.append(LS);
    }
  }

  @Override
  public void close() {
  }

  @Override
  public void endOfSequence() {
    Diagnostic.developerLog("count(l)=" + mEll + " count(1-l)=" + mOneMinusEll);
    Diagnostic.developerLog("count(mu)=" + mMu + " count(1-mu)=" + mOneMinusMu);
    Diagnostic.userLog("hat l=" + ((mEll + 1) / (mEll + mOneMinusEll + 2)));
    Diagnostic.userLog("hat mu=" + ((mMu + 1) / (mMu + mOneMinusMu + 2)));
  }

  @Override
  public boolean integrity() {
    Exam.assertTrue(mQHaploidFactory != null || mQDiploidFactory != null);
    if (mQHaploidFactory != null) {
      checkQ(mQHaploidFactory.somaticQ(mParams.somaticParams().somaticRate()));
    }
    if (mQDiploidFactory != null) {
      checkQ(mQDiploidFactory.somaticQ(mParams.somaticParams().somaticRate()));
    }
    return true;
  }

  private void checkQ(double[][] qa) {
    for (final double[] a : qa) {
      double sum = 0.0;
      for (final double q : a) {
        sum += q;
        Exam.assertTrue(q >= 0.0 && q <= 1.0 && !Double.isNaN(q));
      }
      Exam.assertEquals(1.0, sum, 0.0000001);
    }
  }

}
