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

package com.rtg.variant.bayes.snp;

import java.util.ArrayList;
import java.util.List;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;
import com.rtg.variant.bayes.AlleleBalanceProbability;
import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.Hypotheses;
import com.rtg.variant.bayes.Model;
import com.rtg.variant.bayes.ModelFactory;
import com.rtg.variant.bayes.ModelInterface;

/**
 * @param <H> class of <code>Hypotheses</code>. Needed so <code>ModelCancerFactory</code> can insist on a <code>HypothesesCancer</code>.
 * @param <D> description type
 */
@TestClass("com.rtg.variant.bayes.snp.ModelSnpFactoryTest")
public abstract class ModelCommonFactory<D extends Description, H extends Hypotheses<D>> extends IntegralAbstract implements ModelFactory<D, H> {

  protected H mHypothesisUnknown = null;
  protected final List<H> mHypothesesCache = new ArrayList<>();
  private final AlleleBalanceProbability mAlleleBalance;

  /**
   * @param alleleBalance allele balance probability implementation
   */
  protected ModelCommonFactory(AlleleBalanceProbability alleleBalance) {
    mAlleleBalance = alleleBalance;
  }

  @Override
  public ModelInterface<D> make(final int ref) {
    final Hypotheses<D> hyp = defaultHypotheses(ref);
    return makeModel(hyp);
  }

  @Override
  public H defaultHypotheses(int ref) {
    //return mHypothesesCache.get(ref);
    return ref == Hypotheses.NO_HYPOTHESIS ? mHypothesisUnknown : mHypothesesCache.get(ref);
  }

  protected ModelInterface<D> makeModel(final Hypotheses<D> hyp) {
    return new Model<>(hyp, new StatisticsSnp(hyp.description()), mAlleleBalance);
  }

  @Override
  public boolean globalIntegrity() {
    for (int i = 0; i < mHypothesesCache.size(); ++i) {
      Exam.assertEquals(i, mHypothesesCache.get(i).reference());
    }
    return integrity();
  }

  public AlleleBalanceProbability getAlleleBalance() {
    return mAlleleBalance;
  }

  @Override
  public boolean integrity() {
    return true;
  }

}
