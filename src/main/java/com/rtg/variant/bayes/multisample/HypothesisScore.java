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

import com.rtg.variant.VariantSample;
import com.rtg.variant.bayes.GenotypeMeasure;


/**
 * Holds a single hypothesis along with its scores.
 */
public class HypothesisScore {

  private final GenotypeMeasure mMeasure;
  private VariantSample.DeNovoStatus mDenovo = VariantSample.DeNovoStatus.UNSPECIFIED;
  private Double mDeNovoPosterior = null;

  /**
   * Create a scored hypothesis with no auxiliary score (set to NaN).
   * @param measure the genotype measure for the sample
   */
  public HypothesisScore(final GenotypeMeasure measure) {
    mMeasure = measure;
  }

  /**
   * @return the integer code of an hypothesis (you need an <code>Hypotheses</code> object to interpret it).
   */
  public int hypothesis() {
    return mMeasure.best();
  }

  /**
   * @return the primary score for the hypothesis (normalized in natural log space and of form <code>ln(p/(1-p))</code>).
   */
  public GenotypeMeasure genotypeMeasure() {
    return mMeasure;
  }

  /**
   * @return the score for the best hypothesis in (normalized in natural log space and of form <code>ln(p/(1-p))</code>)
   */
  public double posterior() {
    return mMeasure.arithmetic().poss2Ln(mMeasure.bestPosterior());
  }

  /**
   * @return the auxiliary score for the hypothesis (eg this might be the non-identity posterior).
   */
  public double nonIdentityPosterior() {
    return mMeasure.arithmetic().poss2Ln(mMeasure.nonIdentityPosterior());
  }

  /**
   * The call is de novo
   * @param value set value for de novo flag
   */
  public void setDenovo(VariantSample.DeNovoStatus value) {
    mDenovo = value;
  }

  /**
   * Assign a de novo score
   * @param score probability the best explanation is a de novo mutation (in log space of form <code>ln(p/(1-p)</code>)
   */
  public void setDeNovoPosterior(double score) {
    mDeNovoPosterior = score;
  }

  /**
   * @return get if call is de novo
   */
  public VariantSample.DeNovoStatus isDeNovo() {
    return mDenovo;
  }

  /**
   * @return the score indicating the chance this is a de novo mutation
   */
  public Double getDeNovoPosterior() {
    return mDeNovoPosterior;
  }
}
