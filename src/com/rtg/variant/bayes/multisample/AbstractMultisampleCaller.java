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

package com.rtg.variant.bayes.multisample;

import java.io.IOException;
import java.util.List;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.mode.DnaUtils;
import com.rtg.variant.Variant;
import com.rtg.variant.VariantLocus;
import com.rtg.variant.VariantOutputLevel;
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

    final boolean onlyRefCoverage = Utils.hasOnlyRefCoverage(models);
    if (getParams().callLevel() != VariantOutputLevel.ALL && onlyRefCoverage) {
      return null;
    }

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
    if (hypotheses.haploid().reference() == -1) {
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
   * @param <T> the hypothesis typea
   * @return a <code>ComparisonResult</code> containing the evaluated samples
   */
  protected abstract <D extends Description, T extends HypothesesPrior<D>> ComparisonResult makeSamples(List<ModelInterface<?>> models, HaploidDiploidHypotheses<T> hypotheses);


  protected abstract VariantParams getParams();

  @Override
  public void close() throws IOException {
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
