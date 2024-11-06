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

package com.rtg.variant.bayes.multisample;

import java.util.List;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.mode.DnaUtils;
import com.rtg.variant.Variant;
import com.rtg.variant.VariantLocus;
import com.rtg.variant.VariantOutputOptions;
import com.rtg.variant.VariantParams;
import com.rtg.variant.VariantSample;
import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.Hypotheses;
import com.rtg.variant.bayes.ModelInterface;
import com.rtg.variant.bayes.snp.HypothesesPrior;
import com.rtg.variant.util.VariantUtils;

/**
 */
@TestClass({"com.rtg.variant.bayes.multisample.family.FamilyCallerTest", "com.rtg.variant.bayes.multisample.population.PopulationCallerTest"})
public abstract class AbstractMultisampleCaller implements MultisampleJointCaller {

  /**
   * Create a sample
   *
   * @param hypotheses the hypotheses the sample was generated from
   * @param hypScore the score generated for the model
   * @param b the model
   * @param params parameters governing output
   * @return a newly initialised sample
   */
  public static VariantSample createSample(Hypotheses<?> hypotheses, HypothesisScore hypScore, ModelInterface<?> b, VariantOutputOptions params) {
    final int hyp = hypScore.hypothesis();
    //System.err.println(hypScore.score());
    final boolean isReference = hypotheses.reference() == hyp;

    final VariantSample sample = new VariantSample(hypotheses.ploidy(), hypotheses.name(hyp), isReference, hypScore.genotypeMeasure(), hypScore.isDeNovo(), hypScore.getDeNovoPosterior());

    b.statistics().addCountsToSample(sample, b, params);
    return sample;
  }

  protected void addFilters(String templateName, VariantParams params, List<ModelInterface<?>> models, Variant v) {
    final boolean overCoverage = params.maxCoverageFilter() != null && Utils.maxCoverage(models) > params.maxCoverageFilter().thresholdSingle(templateName);
    if (overCoverage) {
      v.addFilter(Variant.VariantFilter.COVERAGE);
    } else if (params.maxAmbiguity() != null && Utils.meanAmbiguityRatio(models) > params.maxAmbiguity()) {
      v.addFilter(Variant.VariantFilter.AMBIGUITY);
    }
  }

  @Override
  public <D extends Description, T extends HypothesesPrior<D>> Variant makeCall(String templateName, int position, int endPosition, byte[] ref, List<ModelInterface<?>> models, HaploidDiploidHypotheses<T> hypotheses) {

    final ComparisonResult result = makeSamples(models, hypotheses);
    if (result == null) {
      return null;
    }

    final String refAllele = DnaUtils.bytesToSequenceIncCG(ref, position, endPosition - position);
    final VariantLocus locus = new VariantLocus(templateName, position, endPosition, refAllele, VariantUtils.getPreviousRefNt(ref, position));
    final Variant v = new Variant(locus, result.getSamples());
    if (result.isInteresting()) {
      v.setInteresting();
    }
    if (hypotheses.haploid().reference() == Hypotheses.NO_HYPOTHESIS) {
      if (Utils.totalCoverage(models) < Utils.MIN_DEPTH_FOR_N_CALL) {
        return null;
      }
      v.setInvalidRef();
    } else if (result.hasNonIdentity()) {
      v.setNonIdentityPosterior(result.getNonIdentityPosterior());
    }
    addFilters(templateName, getParams(), models, v);
    return v;
  }

  /**
   * Evaluate the models and produce some calls
   * @param models the models to evaluate
   * @param hypotheses hypotheses containing current priors
   * @param <D> the description
   * @param <T> the hypothesis type
   * @return a <code>ComparisonResult</code> containing the evaluated samples
   */
  protected abstract <D extends Description, T extends HypothesesPrior<D>> ComparisonResult makeSamples(List<ModelInterface<?>> models, HaploidDiploidHypotheses<T> hypotheses);


  protected abstract VariantParams getParams();

  @Override
  public void close() {
  }

  @Override
  public void endOfSequence() {
  }

  /**
   * Class wrapping the result of making calls for a set of models
   */
  public static class ComparisonResult {
    private final VariantSample[] mSamples;
    private final boolean mInteresting;
    private final Double mNonIdentityPosterior;

    /**
     * Constructor
     * @param interesting is the result interesting
     * @param samples evaluated samples
     * @param nonIdentityPosterior if there is a score associated that should be used in VCF <code>QUAL</code> field
     */
    public ComparisonResult(boolean interesting, VariantSample[] samples, Double nonIdentityPosterior) {
      mInteresting = interesting;
      mSamples = samples;
      mNonIdentityPosterior = nonIdentityPosterior;
    }

    /**
     * @return the array of samples being created
     */
    public VariantSample[] getSamples() {
      return mSamples;
    }

    /**
     * @return true if the call is interesting and worth outputting
     */
    public boolean isInteresting() {
      return mInteresting;
    }

    /**
     * @return true if a score was generated
     */
    public boolean hasNonIdentity() {
      return mNonIdentityPosterior != null;
    }

    /**
     * @return the generated score
     */
    public double getNonIdentityPosterior() {
      if (mNonIdentityPosterior == null) {
        throw new IllegalStateException("you should check for a non-identity posterior score before calling this method");
      }
      return mNonIdentityPosterior;
    }

  }
}
