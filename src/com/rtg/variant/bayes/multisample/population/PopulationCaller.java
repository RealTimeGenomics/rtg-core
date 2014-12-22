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

package com.rtg.variant.bayes.multisample.population;


import java.util.List;

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
import com.rtg.variant.bayes.snp.HypothesesPrior;

/**
 * Generates the family comparison. Assumes the output will come in with
 * father and mother bayesians in first two array positions.
 */
public class PopulationCaller extends AbstractMultisampleCaller implements MultisampleJointScorer {

  private final MultisampleJointScorer mFamilyCaller;

  private final VariantParams mParams;

  /**
   * @param params variant params
   */
  public PopulationCaller(VariantParams params) {
    this(params, null);
  }

  /**
   * @param params variant params
   * @param familyCaller optional caller used for nuclear families
   */
  public PopulationCaller(VariantParams params, MultisampleJointScorer familyCaller) {
    mParams = params;
    mFamilyCaller = familyCaller;
  }

  /**
   * Generate best scores for a population, without EM iteration.
   *
   * @param <D> the type of the description.
   * @param <T> the type of the hypotheses prior.
   * @param models input models to call from
   * @param priorContainer container for prior related information
   * @return the scores.
   */
  @Override
  public <D extends Description, T extends HypothesesPrior<D>> HypothesisScores getBestScores(List<ModelInterface<?>> models, PriorContainer<T> priorContainer) {

    /* Math to work out QUAL field value from non-identity posterior of each sample
     * qual = 1 - product(1 - p_i), where p_i is probability for sample i that call is non-identity
     * assumes non-identity-posterior value for each sample is
     * nip = ln(p / (1 - p))
     * => p = e^nip / (e^nip + 1)
     * So qual as a probability in terms of non-identity-posterior is
     * qual = 1 - product(1 - (e^nip_i / (e^nip_i + 1)))
     *      = 1 - product(1 / (e^nip_i + 1))
     *      = 1 - e^sum(ln(1 / (e^nip_i + 1))
     *      = 1 - e^sum(-ln(e^nip_i + 1))
     *      = 1 - e^-s, where s = sum(ln(e^nip_i + 1))
     * Then qual as a posterior is
     * Q = ln(qual / (1 - qual))
     *   = ln((1 - e^-s) / (1 - (1 - e^-s)))
     *   = ln(1/e^-s - 1)
     *   = ln(e^s - 1)
     *
     * For ref equals calls the probability numbers are q = 1 - p,
     * => nip = -ip
     *
     * sample.getNonIdentityPosterior produces +ve numbers for non-ref calls, and -ve numbers for ref= calls (identity posterior)
     */

    // First do family calls, then substitute in singleton calls for any gaps
    final HypothesisScore[] allCalls;
    boolean isInteresting = false;
    double sumNips = 0.0; // sum of non-identity posteriors
    HypothesisScores familyCalls = null;
    if (mFamilyCaller != null) {
      familyCalls = mFamilyCaller.getBestScores(models, priorContainer);
      allCalls = familyCalls.getScores();
      isInteresting = familyCalls.isInteresting();

      final double nonIdentityPosterior = familyCalls.getNonIdentityPosterior();
      if (!Double.isNaN(nonIdentityPosterior)) {
//        System.err.print(" " + nonIdentityPosterior);
        sumNips += MathUtils.logExpPlus1(nonIdentityPosterior); // Math.log(Math.exp(nip) + 1), nip = -ip for ref= calls
      }
//      System.err.println();
    } else {
      allCalls = new HypothesisScore[models.size()];
    }
    for (int i = 0; i < allCalls.length; i++) {
      if (allCalls[i] == null) {
        // Get singleton call
        final ModelInterface<?> model = models.get(i);
        final T hyp = priorContainer.getHypotheses().get(model);
        final HypothesisScore best = model.best(hyp);
        allCalls[i] = best;
        if (best != null) { // May be null for female on Y chromosome
          if (hyp.reference() != best.hypothesis()) {
            isInteresting = true;
          }
          if (!Double.isNaN(best.nonIdentityPosterior())) {
//            System.err.print(" " + best.nonIdentityPosterior());
            sumNips += MathUtils.logExpPlus1(best.nonIdentityPosterior()); // Math.log(Math.exp(nip) + 1), nip = -ip for ref= calls
          }
//          System.err.println();
        }
      }
    }
    return new HypothesisScores(allCalls, isInteresting, MathUtils.logExpMinus1(sumNips), familyCalls != null ? familyCalls.getBs() : null);
  }

  @Override
  protected <D extends Description, T extends HypothesesPrior<D>> ComparisonResult makeSamples(List<ModelInterface<?>> models, HaploidDiploidHypotheses<T> hypotheses) {
    final BContainer[] bs = mFamilyCaller != null ? mFamilyCaller.makeInitialBs(models) : null;
    final HypothesisScores calls;

    if (mParams.maxEmIterations() != 0) {
      calls = new EmAlgorithm(new HwEstimator(this), mParams.maxEmIterations()).getBestScores(models, new PriorContainer<>(hypotheses, bs));
    } else {
      calls = getBestScores(models, new PriorContainer<>(hypotheses, bs));
    }

    if (!calls.isInteresting() && mParams.callLevel() != VariantOutputLevel.ALL) {
      return null;
    }

    int sampleNumber = 0;
    final VariantSample[] samples = new VariantSample[models.size()]; //father, mother, children
    for (final HypothesisScore score : calls.getScores()) {
      if (score != null) { // May be null for female on Y chromosome
        final Hypotheses<?> hyp = models.get(sampleNumber).hypotheses();
        final VariantSample sample = createSample(hyp, score, models.get(sampleNumber), mParams);
        samples[sampleNumber] = sample;
      }
      sampleNumber++;
    }

    return new ComparisonResult(calls.isInteresting(), samples, calls.getNonIdentityPosterior()); // qual will be negative if all samples are ref= calls
  }

  @Override
  public BContainer[] makeInitialBs(List<ModelInterface<?>> models) {
    return null;
  }

  @Override
  protected VariantParams getParams() {
    return mParams;
  }
}
