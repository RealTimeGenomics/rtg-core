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

package com.rtg.variant.bayes;

import com.rtg.variant.bayes.complex.DescriptionComplex;
import com.rtg.variant.bayes.multisample.HypothesisScore;
import com.rtg.variant.bayes.snp.DescriptionNone;
import com.rtg.variant.bayes.snp.HypothesesNone;
import com.rtg.variant.bayes.snp.HypothesesPrior;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;

/**
 * Empty model.
 */
public final class ModelNone<D extends Description> implements ModelInterface<D> {

  /** Model for case where no evidence. */
  public static final ModelNone<Description> SINGLETON = new ModelNone<>(HypothesesNone.SINGLETON);

  /** Model for case where no evidence. */
  public static final ModelNone<DescriptionComplex> SINGLETON_COMPLEX = new ModelNone<>(HypothesesNone.SINGLETON_COMPLEX);

  private final Hypotheses<D> mHypotheses;

  private ModelNone(Hypotheses<D> hyp) {
    mHypotheses = hyp;
  }

  @Override
  public Hypotheses<D> hypotheses() {
    return mHypotheses;
  }

  @Override
  public int size() {
    return mHypotheses.size();
  }

  @Override
  public String name(int i) {
    return mHypotheses.name(i);
  }

  @Override
  public Description description() {
    return null;
  }

  @Override
  public int reference() {
    return mHypotheses.reference();
  }


  @Override
  public boolean haploid() {
    return mHypotheses.haploid();
  }


  @Override
  public void increment(EvidenceInterface evidence) {
    //Occasionally this can happen when overlapping into PAR regions on Y
    //Diagnostic.developerLog("Evidence for sequence that should not exist for this individual: " + evidence);
  }

  @Override
  public void freeze() {
  }

  @Override
  public void statistics(StringBuilder sb, final HypothesesPrior<?> hypotheses) {
  }

  @Override
  public Statistics<?> statistics() {
    return new StatisticsInt(DescriptionNone.SINGLETON);
  }

  @Override
  public double posteriorLn0(int hyp) {
    return 0;
  }

  @Override
  public HypothesisScore best(HypothesesPrior<?> descriptionComplexHypotheses) {
    return null;
  }

  @Override
  public AlleleBalanceProbability alleleBalanceProbability() {
    return null;
  }

  @Override
  public PossibilityArithmetic arithmetic() {
    return mHypotheses.arithmetic();
  }

  @Override
  public double p(int code) {
    return arithmetic().one();
  }

  @Override
  public boolean isNormalized() {
    return true;
  }

  @Override
  public Factor<D> normalize() {
    return this;
  }

  @Override
  public ModelNone<D> copy() {
    return new ModelNone<>(mHypotheses);
  }
}
