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


import java.io.IOException;
import java.util.List;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.reference.Ploidy;
import com.rtg.util.Utils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;
import com.rtg.variant.Variant;
import com.rtg.variant.Variant.VariantFilter;
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
import com.rtg.variant.bayes.Statistics;
import com.rtg.variant.bayes.multisample.HaploidDiploidHypotheses;
import com.rtg.variant.bayes.multisample.MultisampleJointCaller;
import com.rtg.variant.bayes.snp.HypothesesPrior;
import com.rtg.variant.bayes.snp.StatisticsSnp;
import com.rtg.variant.dna.DNARange;
import com.rtg.variant.util.VariantUtils;

/**
 */
@TestClass(value = {"com.rtg.variant.bayes.multisample.cancer.ContaminatedSomaticCallerTest", "com.rtg.variant.bayes.multisample.cancer.PureSomaticCallerTest"})
public abstract class AbstractSomaticCaller extends IntegralAbstract implements MultisampleJointCaller {

  protected final double[][] mQHaploid;
  protected final double[][] mQDiploid;
  private final VariantParams mParams;

  /**
   * @param qHaploid haploid numbers
   * @param qDiploid diploid numbers
   * @param params variant params
   */
  public AbstractSomaticCaller(final double[][] qHaploid, double[][] qDiploid, VariantParams params) {
    mQHaploid = qHaploid;
    mQDiploid = qDiploid;
    mParams = params;
  }

  /**
   * Construct an appropriate posterior. Differs in the contaminated and non-contaminated case.
   * @param normal bayesian for the normal genome.
   * @param cancer bayesian for the cancer genome.
   * @param hypotheses the hypotheses containing priors.
   * @return the posterior.
   */
  protected abstract AbstractPosterior makePosterior(final ModelInterface<?> normal, final ModelInterface<?> cancer, HypothesesPrior<?> hypotheses);

  private VariantSample setCallValues(GenotypeMeasure posterior, int cat, Hypotheses<?> hypotheses, ModelInterface<?> model, VariantOutputOptions params, Ploidy ploidy) {
    final VariantSample sample = new VariantSample(ploidy, hypotheses.name(cat), hypotheses.reference() == cat, posterior, VariantSample.DeNovoStatus.UNSPECIFIED, null);
    model.statistics().addCountsToSample(sample, model, params);
    return sample;
  }

  private double loh(final Hypotheses<?> hypotheses, final int normal, final int cancer) {
    final Code code = hypotheses.code();
    if (!hypotheses.code().homozygous(normal) && code.homozygous(cancer)) {
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
  private double mAlpha = 0;
  private double mOneMinusAlpha = 0;

  private void updateParameterEstimationCounts(final Ploidy normalPloidy, final Code code, final int normal, final int cancer, final double loh) {
    // See multiscoring.tex for an explanation.
    // A more thorough treatment should probably increment mEll and mOneMinusEll in proportion to loh
    // but we often have loh of 0 when we haven't actually made a call, so this could be tricky.
    // The theory really only covers the SNP cases, but here we just "make it work" for complex
    // cases as well.
    if (loh <= 0) {
      mOneMinusEll++;
      if (normalPloidy == Ploidy.HAPLOID) {
        if (normal == cancer) {
          mOneMinusMu++;
        } else {
          mMu++;
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
      mEll++;
      if (code.a(normal) == cancer || code.bc(normal) == cancer) {
        mOneMinusMu++;
      } else {
        mMu++;
      }
    }
  }

  private void updateContaminationCounts(final Ploidy normalPloidy, final Code code, ModelInterface<?> modelCancer, final int normal, final int cancer, final double loh) {
    // See multiscoring.tex for an explanation.
    // This doesn't examine all possible cases due to difficulty of integration
    final Statistics<?> statsCancer = modelCancer.statistics();
    if (!(statsCancer instanceof StatisticsSnp)) {
      return;
    }
    final AlleleStatistics<?> counts = ((StatisticsSnp) statsCancer).counts();
    if (loh <= 0) {
      if (normalPloidy == Ploidy.HAPLOID) {
        if (normal != cancer) {
          mAlpha += counts.count(normal);
          mOneMinusAlpha += counts.count(cancer);
        }
      } else {
        final int na = code.a(normal);
        final int nb = code.bc(normal);
        final int ca = code.a(cancer);
        final int cb = code.bc(cancer);
        if (na == nb) {
          if (ca != na && cb != na) {
            mAlpha += counts.count(na);
            mOneMinusAlpha += 0.5 * counts.count(ca);
            mOneMinusAlpha += 0.5 * counts.count(cb);
          }
        } else if (ca == cb) {
          if (na != ca && nb != ca) {
            mAlpha += 0.5 * counts.count(na);
            mAlpha += 0.5 * counts.count(nb);
            mOneMinusAlpha += counts.count(ca);
          }
        } else {
          // na != nb && ca != cb;
          if (ca != na && ca != nb && cb != na && cb != nb) {
            mAlpha += 0.5 * counts.count(na);
            mAlpha += 0.5 * counts.count(nb);
            mOneMinusAlpha += 0.5 * counts.count(ca);
            mOneMinusAlpha += 0.5 * counts.count(cb);
          }
        }
      }
    }
  }

  @Override
  public <D extends Description, T extends HypothesesPrior<D>> Variant makeCall(String templateName, int position, int endPosition, byte[] ref, List<ModelInterface<?>> models, HaploidDiploidHypotheses<T> normalHypotheses) {
    assert DNARange.valid(ref[position], DNARange.N, DNARange.T);
    assert models.size() == 2;

    // TODO: with ALL output we should make a call here
    if (ref[position] == DNARange.N) {
      return null;
    }
    final ModelInterface<?> modelNormal = models.get(0);
    final ModelInterface<?> modelCancer = models.get(1);
    final HypothesesPrior<?> hypotheses = normalHypotheses.get(modelNormal);
    final Code code = hypotheses.code();

    // Boring keeps track of any call which is not-interesting to the cancer mode.
    // Such calls must never have isInteresting set to true when they are returned
    // so that later complex calling does not muck with them.

    boolean boring = modelNormal.statistics().nonRefCount() == 0 && modelCancer.statistics().nonRefCount() == 0;
    if (boring && mParams.callLevel() != VariantOutputLevel.ALL) {
      // Short-circuit case where all evidence matches the reference.
      updateParameterEstimationCounts(null, code, 0, 0, 0);
      return null;
    }

    final AbstractPosterior posterior = makePosterior(modelNormal, modelCancer, hypotheses);
    final boolean sameCall = posterior.isSameCall();
    final String cause;
    final int bestNormal = posterior.bestNormal();
    final int bestCancer = posterior.bestCancer();
    // Simple LOH test based on ploidy of results alone, could be done with Bayesian calculation later
    final double loh = loh(hypotheses, bestNormal, bestCancer);
    // create individual genome calls
    final Ploidy normalPloidy = hypotheses.haploid() ? Ploidy.HAPLOID : Ploidy.DIPLOID;
    final Ploidy cancerPloidy;
    final Hypotheses<?> cancerHyp;
    final boolean doLoh = mParams.lohPrior() > 0;
    if (loh > 0) {
      cancerPloidy = Ploidy.HAPLOID;
      cancerHyp = normalHypotheses.haploid();
    } else {
      cancerPloidy = normalPloidy;
      cancerHyp = hypotheses;
    }
    if (sameCall || (!doLoh && loh > 0)) {
      // This call is either
      //  - the same for normal and cancer samples, but may be different to the reference.
      //  - looks like an LOH event, even though the LOH prior was 0
      // It is not interesting and should not normally be output,
      if (mParams.callLevel() != VariantOutputLevel.ALL) {
        return null;
      }
      boring = true;
      cause = "NONE";
    } else {
      //TODO cannot be cause if it is equal to reference
      cause = cancerHyp.name(bestCancer);
    }

    // set the combined call values
    final double ratio = posterior.posteriorScore();
    if (ratio < mParams.threshold()) {
      return null;
    }

    final VariantSample normalSample = setCallValues(posterior.normalMeasure(), bestNormal, hypotheses, modelNormal, mParams, normalPloidy);
    final VariantSample cancerSample = setCallValues(posterior.cancerMeasure(), bestCancer, cancerHyp, modelCancer, mParams, cancerPloidy);

    final String refAllele = hypotheses.description().name(hypotheses.reference());
    final VariantLocus locus = new VariantLocus(templateName, position, endPosition, refAllele, VariantUtils.getPreviousRefNt(ref, position));

    final Variant v = new Variant(locus, sameCall, normalSample, cancerSample);
    v.setNormalCancerScore(posterior.ncScore());

    if (modelNormal.statistics().overCoverage(mParams, templateName) || modelCancer.statistics().overCoverage(mParams, templateName)) {
      v.addFilter(VariantFilter.COVERAGE);
      boring = true;
    } else if (modelNormal.statistics().ambiguous(mParams) || modelCancer.statistics().ambiguous(mParams)) {
      v.addFilter(VariantFilter.AMBIGUITY);
    }
    v.setPossibleCause(cause);
    v.setPossibleCauseScore(ratio);
    if (doLoh) {
      v.setLoh(loh);
    }
    if (!boring) {
      v.setInteresting();
    }
    updateParameterEstimationCounts(normalPloidy, code, bestNormal, bestCancer, loh);
    if (!v.isFiltered()) {
      updateContaminationCounts(normalPloidy, code, modelCancer, bestNormal, bestCancer, loh);
    }
    return v;
  }

  @Override
  public void toString(StringBuilder sb) {
    final double[][] qa = mQHaploid != null ? mQHaploid : mQDiploid;
    sb.append("length=").append(qa.length).append(LS);
    for (final double[] aQa : qa) {
      for (final double anAQa : aQa) {
        sb.append(Utils.realFormat(anAQa, 3)).append(" ");
      }
      sb.append(LS);
    }
  }

  @Override
  public void close() throws IOException {
  }

  @Override
  public void endOfSequence() {
    Diagnostic.developerLog("count(l)=" + mEll + " count(1-l)=" + mOneMinusEll);
    Diagnostic.developerLog("count(mu)=" + mMu + " count(1-mu)=" + mOneMinusMu);
    Diagnostic.developerLog("count(alpha)=" + mAlpha + " count(1-alpha)=" + mOneMinusAlpha);
    Diagnostic.userLog("hat l=" + ((mEll + 1) / (mEll + mOneMinusEll + 2)));
    Diagnostic.userLog("hat mu=" + ((mMu + 1) / (mMu + mOneMinusMu + 2)));
    Diagnostic.userLog("hat alpha=" + ((mAlpha + 1) / (mAlpha + mOneMinusAlpha + 2)));
  }

  @Override
  public boolean integrity() {
    Exam.assertTrue(mQHaploid != null || mQDiploid != null);
    if (mQHaploid != null) {
      checkQ(mQHaploid);
    }
    if (mQDiploid != null) {
      checkQ(mQDiploid);
    }
    return true;
  }

  private void checkQ(double[][] qa) {
    final int length = qa.length;
    for (final double[] aQa : qa) {
      Exam.assertEquals(length, aQa.length);
      double sum = 0.0;
      for (final double q : aQa) {
        sum += q;
        Exam.assertTrue(q >= 0.0 && q <= 1.0 && !Double.isNaN(q));
      }
      Exam.assertEquals(1.0, sum, 0.0000001);
    }
  }

}
