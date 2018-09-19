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
package com.rtg.variant.bayes.multisample.family;


import java.util.List;

import com.rtg.relation.Family;
import com.rtg.variant.Variant;
import com.rtg.variant.Variant.VariantFilter;
import com.rtg.variant.VariantLocus;
import com.rtg.variant.VariantOutputLevel;
import com.rtg.variant.VariantParams;
import com.rtg.variant.VariantSample;
import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.Hypotheses;
import com.rtg.variant.bayes.ModelInterface;
import com.rtg.variant.bayes.multisample.HaploidDiploidHypotheses;
import com.rtg.variant.bayes.multisample.HypothesisScore;
import com.rtg.variant.bayes.multisample.MultisampleJointCaller;
import com.rtg.variant.bayes.multisample.Utils;
import com.rtg.variant.bayes.snp.HypothesesPrior;
import com.rtg.mode.DNARange;
import com.rtg.variant.util.VariantUtils;

/**
 * Generates the diseased family comparison. Assumes the models will come in with
 * father and mother in first two array positions.
 */
public class DiseasedFamilyCaller implements MultisampleJointCaller {

  private static final boolean USE_SLOW_IMPL = false; //Boolean.valueOf(System.getProperty("rtg.family.disease.use-slow", "false"));
  private static final double POSTERIOR_THRESHOLD_NOT_INTERESTING = 5.0;

  private final VariantParams mParams;
  private final Family mFamily;
  private final double mNoDiseasePrior;

  /**
   * @param params variant params
   * @param family the family (including who has the disease)
   * @param noDiseasePrior prior that there is no disease explanation at this position
   */
  //TODO add flag for recessive vs dominant
  public DiseasedFamilyCaller(VariantParams params, Family family, double noDiseasePrior) {
    mParams = params;
    mFamily = family;
    mNoDiseasePrior = noDiseasePrior;
  }

  static boolean shortCircuit(final List<ModelInterface<?>> models, HaploidDiploidHypotheses<?> hypotheses) {
    // If everyone agrees on best hypothesis there cannot be a disease explanation at this position.
    // In practice we have to allow for lower weight hypothesis coming up when considered
    // in combination, hence we only short-circuit for large enough posteriors as well.
    final HypothesisScore fatherBest = models.get(0).best(hypotheses.get(models.get(0)));
    if (fatherBest.posterior() > POSTERIOR_THRESHOLD_NOT_INTERESTING) {
      final int fName = fatherBest.hypothesis();
      boolean ok = true;
      for (int k = 1; k < models.size(); ++k) {
        final HypothesisScore bc = models.get(k).best(hypotheses.get(models.get(k)));
        if (bc.hypothesis() != fName || bc.posterior() <= POSTERIOR_THRESHOLD_NOT_INTERESTING) {
          ok = false;
          break;
        }
      }
      if (ok) {
        return true;
      }
    }
    return false;
  }

  @Override
  public <D extends Description, T extends HypothesesPrior<D>> Variant makeCall(String templateName, int position, int endPosition, byte[] ref, List<ModelInterface<?>> models, HaploidDiploidHypotheses<T> hypotheses) {
    if (ref[position] == DNARange.N) {
      return null;
    }
    if (!USE_SLOW_IMPL && mParams.callLevel() != VariantOutputLevel.ALL && shortCircuit(models, hypotheses)) {
      return null;
    }

    // We assume the bayesians will have father and mother first
    final ModelInterface<?> fb = models.get(Family.FATHER_INDEX);
    final ModelInterface<?> mb = models.get(Family.MOTHER_INDEX);
    final boolean overCoverage = mParams.maxCoverageFilter() != null && Utils.maxCoverage(models) > mParams.maxCoverageFilter().thresholdSingle(templateName);
    final Hypotheses<?> commonHypotheses = hypotheses.get(models.get(0)); // XXX this doesn't appear to support sex
    final byte refNt = (byte) (ref[position] - 1);
    final HypothesesDisease diseaseHypotheses = new HypothesesDisease(commonHypotheses.description(), mNoDiseasePrior, refNt);
    //TODO use flag for dominant vs recessive.
    final DiseasedFamilyPosterior fp = USE_SLOW_IMPL ? new DiseasedFamilyPosterior(mParams.genomePriors(), mFamily, commonHypotheses, diseaseHypotheses, models) : new FastDiseasedFamilyPosterior(mParams.genomePriors(), mFamily, commonHypotheses, diseaseHypotheses, models);
    fp.process();
    final HypothesisScore bd = fp.bestDisease();
    final HypothesisScore bf = fp.bestFather();
    final HypothesisScore bm = fp.bestMother();

    final double ratio = bd.posterior();
    if (ratio < 0) {
      return null;
    }

    if (mParams.callLevel() == VariantOutputLevel.ALL || fp.isInteresting()) {
      // TODO: These outputs are supposed to be ratios but current are not?
      final VariantSample[] samples = new VariantSample[models.size()];
      samples[Family.FATHER_INDEX] = FamilyCaller.createSample(commonHypotheses, bf, fb, mParams);
      samples[Family.MOTHER_INDEX] = FamilyCaller.createSample(commonHypotheses, bm, mb, mParams);

      for (int i = 0; i < models.size() - Family.FIRST_CHILD_INDEX; ++i) {
        final HypothesisScore child = fp.bestChild(i);
        samples[i + Family.FIRST_CHILD_INDEX] = FamilyCaller.createSample(commonHypotheses, child, models.get(i + Family.FIRST_CHILD_INDEX), mParams);
      }
      final String refAllele = commonHypotheses.description().name(commonHypotheses.reference()); // XXX Doesn't support reference == Hypothesis.NO_HYPOTHESIS
      final VariantLocus locus = new VariantLocus(templateName, position, endPosition, refAllele, VariantUtils.getPreviousRefNt(ref, position));
      final Variant v = new Variant(locus, samples);
      v.setDiseasePresenceScore(fp.anyDiseasePosteriorRatio());
      //v.setNonIdentityPosterior(fp.anyDiseasePosteriorRatio());
      if (overCoverage) {
        v.addFilter(VariantFilter.COVERAGE);
      } else if (mParams.maxAmbiguity() != null && Utils.meanAmbiguityRatio(models) > mParams.maxAmbiguity()) {
        v.addFilter(VariantFilter.AMBIGUITY);
      }
      v.setPossibleCause(diseaseHypotheses.name(bd.hypothesis()));
      v.setPossibleCauseScore(ratio);
      return v;  //previously only returned 'true' if was interesting. This method however does its own thing for ALL mode. Possible difference in operation due to this
    } else {
      return null;
    }
  }

  @Override
  public void close() {
  }

  @Override
  public void endOfSequence() {
  }

}
