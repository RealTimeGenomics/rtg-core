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
    for (int i = 0; i < mHypothesesCache.size(); i++) {
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
