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

package com.rtg.variant.bayes.multisample.singleton;


import java.util.List;

import com.rtg.launcher.GlobalFlags;
import com.rtg.mode.DnaUtils;
import com.rtg.reference.Ploidy;
import com.rtg.util.MathUtils;
import com.rtg.variant.Variant;
import com.rtg.variant.VariantLocus;
import com.rtg.variant.VariantOutputLevel;
import com.rtg.variant.VariantParams;
import com.rtg.variant.VariantSample;
import com.rtg.variant.bayes.AlleleStatistics;
import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.ModelInterface;
import com.rtg.variant.bayes.multisample.HaploidDiploidHypotheses;
import com.rtg.variant.bayes.multisample.HypothesisScore;
import com.rtg.variant.bayes.multisample.MultisampleJointCaller;
import com.rtg.variant.bayes.multisample.Utils;
import com.rtg.variant.bayes.snp.HypothesesPrior;
import com.rtg.variant.util.VariantUtils;

/**
 */
public class SingletonCaller implements MultisampleJointCaller {

  private static final int MIN_VAC = GlobalFlags.getIntegerValue(GlobalFlags.VARIANT_MIN_AC);
  private static final double MIN_VAF = GlobalFlags.getDoubleValue(GlobalFlags.VARIANT_MIN_AF);

  private final VariantParams mParams;
  private final double mInterestingThreshold;

  /**
   * @param params variant params
   */
  public SingletonCaller(VariantParams params) {
    mParams = params;
    mInterestingThreshold = mParams.interestingThreshold() * MathUtils.LOG_10;
  }

  @Override
  public <D extends Description, T extends HypothesesPrior<D>> Variant makeCall(String templateName, int position, int endPosition, byte[] ref, List<ModelInterface<?>> models, HaploidDiploidHypotheses<T> hypotheses) {
    assert models.size() == 1;

    final boolean onlyRefCoverage = Utils.hasOnlyRefCoverage(models);
    if (mParams.callLevel() != VariantOutputLevel.ALL && onlyRefCoverage) {
      return null;
    }

    final ModelInterface<?> model = models.get(0);
    final T hyp = hypotheses.get(model);

    final String refAllele = DnaUtils.bytesToSequenceIncCG(ref, position, endPosition - position);
    final VariantLocus locus = new VariantLocus(templateName, position, endPosition, refAllele, VariantUtils.getPreviousRefNt(ref, position));
    final int size = model.size();
    assert size >= 1;

    final AlleleStatistics<?> ac = model.statistics().counts();
    final Description description = ac.getDescription();
    final int refCode = description.indexOf(refAllele);
    double tot = 0;
    double vac = 0;
    int va = -1;
    boolean tied = false;
    for (int i = 0; i < description.size(); i++) {
      final double count = ac.count(i) - ac.error(i);
      tot += count;
      if (i != refCode && count >= vac) {
        tied = count == vac;
        vac = count;
        va = i;
      }
    }
    double vaf = 0;
    if (va != -1) {
      vaf = vac / tot;
    }


    final HypothesisScore best = model.best(hyp);
    //final boolean changed = hyp.reference() != best.hypothesis();

    final boolean interesting = vac > MIN_VAC && vaf > MIN_VAF;
    if (interesting && tied) {
      System.err.println("Two non-ref alleles met VAC/VAF threshold at " + templateName + ":" + (position + 1));
      System.err.println("va = " + va + " (" + description.name(va) + ") vac = " + vac + " vaf = " + vaf + " interesting = " + interesting);
    }


    //final boolean interesting = changed || best.posterior() < mInterestingThreshold;

    if (!interesting && mParams.callLevel() != VariantOutputLevel.ALL) {
      return null;
    }

    final String altName = hyp.name(best.hypothesis());
    final VariantSample sample = new VariantSample(hyp.haploid() ? Ploidy.HAPLOID : Ploidy.DIPLOID, altName, best.hypothesis() == hyp.reference(), best.genotypeMeasure(), VariantSample.DeNovoStatus.UNSPECIFIED, null);
    final Variant v = new Variant(locus, sample);
    if (interesting) {
      v.setInteresting();
      sample.setVariantAllele(description.name(va));
    }
    if (hyp.reference() == -1) {
      if (Utils.totalCoverage(models) < Utils.MIN_DEPTH_FOR_N_CALL) {
        return null;
      }
      v.setInvalidRef();
    } else if (mParams.nonidentityPosterior()) {
      v.setNonIdentityPosterior(best.nonIdentityPosterior());
    }
    model.statistics().addFiltersToVariant(v, templateName, mParams);
    model.statistics().addCountsToSample(sample, model, mParams);

    return v;
  }

  @Override
  public void close() {
  }

  @Override
  public void endOfSequence() {
  }

}
