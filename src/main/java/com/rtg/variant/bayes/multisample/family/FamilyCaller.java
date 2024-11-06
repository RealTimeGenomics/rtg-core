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
package com.rtg.variant.bayes.multisample.family;

import java.util.List;

import com.rtg.launcher.globals.GlobalFlags;
import com.rtg.launcher.globals.CoreGlobalFlags;
import com.rtg.relation.Family;
import com.rtg.util.MathUtils;
import com.rtg.variant.VariantOutputLevel;
import com.rtg.variant.VariantParams;
import com.rtg.variant.VariantSample;
import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.Hypotheses;
import com.rtg.variant.bayes.ModelInterface;
import com.rtg.variant.bayes.multisample.AbstractMultisampleCaller;
import com.rtg.variant.bayes.multisample.HaploidDiploidHypotheses;
import com.rtg.variant.bayes.multisample.HypothesisScore;
import com.rtg.variant.bayes.multisample.HypothesisScores;
import com.rtg.variant.bayes.multisample.MultisampleJointScorer;
import com.rtg.variant.bayes.multisample.PriorContainer;
import com.rtg.variant.bayes.multisample.forwardbackward.BContainer;
import com.rtg.variant.bayes.multisample.forwardbackward.FamilyCallerFB;
import com.rtg.variant.bayes.snp.HypothesesPrior;

/**
 * Generates the family comparison. Assumes the output will come in with
 * father and mother Bayesians in first two array positions.
 */
public class FamilyCaller extends AbstractMultisampleCaller implements MultisampleJointScorer {

  private static final boolean USE_SLOW_IMPL = false; //Boolean.valueOf(System.getProperty("rtg.family.use-slow", "false"));
  // If true, resolve disagreeing calls using the forward backward caller (currently SLOW)
  private static final boolean USE_FB_FALLBACK = GlobalFlags.getBooleanValue(CoreGlobalFlags.FAMILY_CALLER_FALLBACK_FLAG);

  private final VariantParams mParams;
  private final Family[] mFamilies;
  private final FamilyCallerFB mFallbackCaller;

  /**
   * @param params genome priors
   * @param families the family structures
   */
  public FamilyCaller(VariantParams params, Family... families) {
    mParams = params;
    mFamilies = families;
    mFallbackCaller = new FamilyCallerFB(params, families);
  }

  // Put a new score into the current scores array, dealing with cases where we have already called the individual
  static boolean setScore(HypothesisScore[] scores, int index, HypothesisScore newScore) {
    if (scores[index] == null) {
      scores[index] = newScore;
      return true;
    } else {
      // A call has already been placed in for this individual, and due to how families are ordered from
      // highest to lowest in the pedigree, this means that any de novo calls will be placed in first.

      if (scores[index].hypothesis() == newScore.hypothesis()) {
        // If there is agreement in the call, keep the higher set of scores (but make sure we preserve the de novo status).
        if (scores[index].posterior() < newScore.posterior()) {
          if (scores[index].isDeNovo() != VariantSample.DeNovoStatus.UNSPECIFIED) {
            // E.g. germ line de novo. Copy over de novo status so that we do not lose it
            newScore.setDenovo(scores[index].isDeNovo());
            newScore.setDeNovoPosterior(scores[index].getDeNovoPosterior());
          }
          scores[index] = newScore;
        }
        return true;
      } else {
        // These disagree. E.g. the upper call may be a somatic de novo.
        // We'll recall with the FB caller
        return false;
      }
    }
  }

  /**
   * Generate best scores for a family
   *
   * @param <D> the type of the description.
   * @param <T> the type of the hypotheses prior.
   * @param models input models to call from
   * @param priorContainer source of all things Hypothesis and <code>B</code>
   * @return the scores.
   */
  @Override
  public <D extends Description, T extends HypothesesPrior<D>> HypothesisScores getBestScores(List<ModelInterface<?>> models, PriorContainer<T> priorContainer) {
    boolean isInteresting = false;
    double sumNips = 0;
    final HypothesisScore[] scores = new HypothesisScore[models.size()];

    boolean agreeing = true;
    for (final Family family : mFamilies) {
      final AbstractFamilyPosterior fp = USE_SLOW_IMPL ? new FamilyPosterior(family, mParams.genomePriors(), models, priorContainer.getHypotheses()) : new FastFamilyPosterior(family, mParams.genomePriors(), models, priorContainer.getHypotheses());
      final int[] ids = family.getSampleIds();
      agreeing = setScore(scores, ids[Family.FATHER_INDEX], fp.bestFather());
      agreeing &= setScore(scores, ids[Family.MOTHER_INDEX], fp.bestMother());
      for (int i = 0; i < family.numChildren(); ++i) {
        agreeing &= setScore(scores, ids[Family.FIRST_CHILD_INDEX + i], fp.bestChild(i));
      }
      if (USE_FB_FALLBACK && !agreeing) {
        break;
      }
      if (fp.isInteresting()) {
        isInteresting = true;
      }
      final double nip = fp.getNonIdentityPosterior();
      if (!Double.isNaN(nip)) {
        sumNips += MathUtils.logExpPlus1(nip);
      }
    }

    if (USE_FB_FALLBACK && !agreeing) {
      return mFallbackCaller.getBestScores(models, priorContainer);
    }
    final double qual = MathUtils.logExpMinus1(sumNips);
    return new HypothesisScores(scores, isInteresting, qual, priorContainer.getBs());
  }

  @Override
  public BContainer[] makeInitialBs(List<ModelInterface<?>> models) {
    return mFallbackCaller.makeInitialBs(models);
  }

  @Override
  protected <D extends Description, T extends HypothesesPrior<D>> ComparisonResult makeSamples(List<ModelInterface<?>> models, HaploidDiploidHypotheses<T> hypotheses) {
    final HypothesisScores scores = getBestScores(models, new PriorContainer<>(hypotheses, null));

    if (!scores.isInteresting() && mParams.callLevel() != VariantOutputLevel.ALL) {
      return null;
    }

    final VariantSample[] samples = new VariantSample[models.size()];
    for (int i = 0; i < samples.length; ++i) {
      final HypothesisScore score = scores.getScores()[i];
      if (score != null && score.hypothesis() != -1) {
        final Hypotheses<?> childHypotheses = hypotheses.get(models.get(i));
        samples[i] = createSample(childHypotheses, score, models.get(i), mParams);
      }
    }

    return new ComparisonResult(scores.isInteresting(), samples, scores.getNonIdentityPosterior());
  }

  @Override
  protected VariantParams getParams() {
    return mParams;
  }
}
