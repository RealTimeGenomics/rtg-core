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
package com.rtg.variant.bayes.multisample.lineage;

import java.util.Collections;
import java.util.Map;
import java.util.Set;

/**
 * Bridge from model interface factor to factor.
 *
 */
public class ModelFactor extends AbstractFactor implements ToDefaultFactor {

  private final com.rtg.variant.bayes.Factor<?> mUnderlyingFactor;
  private final Variable mVar;

  /**
   * Factor from Bayes model.
   * @param var variable
   * @param model underlying model
   */
  public ModelFactor(final Variable var, final com.rtg.variant.bayes.Factor<?> model) {
    super(model.arithmetic());
    mVar = var;
    mUnderlyingFactor = model;
  }

  @Override
  public double p(Map<Variable, Integer> values) {
    checkScope(values);
    return mUnderlyingFactor.p(values.get(mVar));
  }

  @Override
  public Set<Variable> scope() {
    return Collections.singleton(mVar);
  }

  @Override
  public DefaultFactor asDefault() {
    final double[] poss = new double[mVar.size()];
    for (int i = 0; i < mVar.size(); i++) {
      poss[i] = mUnderlyingFactor.p(i);
    }
    return new DefaultFactor(mUnderlyingFactor.arithmetic(), Collections.singletonList(mVar), poss);
  }
}
