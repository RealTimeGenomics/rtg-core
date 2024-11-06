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
    for (int i = 0; i < mVar.size(); ++i) {
      poss[i] = mUnderlyingFactor.p(i);
    }
    return new DefaultFactor(mUnderlyingFactor.arithmetic(), Collections.singletonList(mVar), poss);
  }
}
