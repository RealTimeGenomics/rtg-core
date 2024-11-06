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

package com.rtg.variant.bayes.multisample.singleton;


import java.util.List;

import com.rtg.mode.DnaUtils;
import com.rtg.reference.Ploidy;
import com.rtg.variant.Variant;
import com.rtg.variant.VariantLocus;
import com.rtg.variant.VariantOutputLevel;
import com.rtg.variant.VariantParams;
import com.rtg.variant.VariantSample;
import com.rtg.variant.bayes.AlleleStatistics;
import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.Hypotheses;
import com.rtg.variant.bayes.ModelInterface;
import com.rtg.variant.bayes.ModelNone;
import com.rtg.variant.bayes.multisample.HaploidDiploidHypotheses;
import com.rtg.variant.bayes.multisample.HypothesisScore;
import com.rtg.variant.bayes.multisample.MultisampleJointCaller;
import com.rtg.variant.bayes.multisample.Utils;
import com.rtg.variant.bayes.multisample.VariantAlleleTrigger;
import com.rtg.variant.bayes.snp.HypothesesPrior;
import com.rtg.variant.util.VariantUtils;

/**
 */
public class SingletonCaller implements MultisampleJointCaller {

  private final VariantParams mParams;
  private final double mInterestingThreshold;
  private final VariantAlleleTrigger mVariantAlleleTrigger;

  /**
   * @param params variant params
   */
  public SingletonCaller(VariantParams params) {
    mParams = params;
    mInterestingThreshold = mParams.interestingThreshold();
    mVariantAlleleTrigger = new VariantAlleleTrigger(params.minVariantAllelicDepth(), params.minVariantAllelicFraction());
  }

  @Override
  public <D extends Description, T extends HypothesesPrior<D>> Variant makeCall(String templateName, int position, int endPosition, byte[] ref, List<ModelInterface<?>> models, HaploidDiploidHypotheses<T> hypotheses) {
    assert models.size() == 1;

    final ModelInterface<?> model = models.get(0);
    if (model instanceof ModelNone<?>) {
      return null;
    }

    final T hyp = hypotheses.get(model);

    final String refAllele = DnaUtils.bytesToSequenceIncCG(ref, position, endPosition - position);
    final VariantLocus locus = new VariantLocus(templateName, position, endPosition, refAllele, VariantUtils.getPreviousRefNt(ref, position));
    final int size = model.size();
    assert size >= 1;

    final AlleleStatistics<?> ac = model.statistics().counts();
    final Description description = ac.getDescription();
    final int va = mVariantAlleleTrigger.getVariantAllele(ac, description, refAllele);

    final HypothesisScore best = model.best(hyp);

    final boolean triggersVa = va != -1;

    final boolean changed = hyp.reference() != best.hypothesis();
    final boolean interesting = changed || best.posterior() < mInterestingThreshold;

    if (!(interesting || triggersVa || mParams.callLevel() == VariantOutputLevel.ALL)) {
      return null;
    }

    final String altName = hyp.name(best.hypothesis());
    final VariantSample sample = new VariantSample(hyp.haploid() ? Ploidy.HAPLOID : Ploidy.DIPLOID, altName, best.hypothesis() == hyp.reference(), best.genotypeMeasure(), VariantSample.DeNovoStatus.UNSPECIFIED, null);
    if (triggersVa) {
      sample.setVariantAllele(description.name(va));
    }
    final Variant v = new Variant(locus, sample);
    if (interesting) {
      v.setInteresting();
    }
    if (hyp.reference() == Hypotheses.NO_HYPOTHESIS) {
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
