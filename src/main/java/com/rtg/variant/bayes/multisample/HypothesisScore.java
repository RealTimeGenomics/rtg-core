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
