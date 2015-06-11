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

}
