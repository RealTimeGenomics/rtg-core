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
package com.rtg.variant.bayes.multisample.cancer;

import java.util.ArrayList;
import java.util.List;

import com.rtg.variant.VariantParams;
import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.Hypotheses;
import com.rtg.variant.bayes.ModelInterface;
import com.rtg.variant.bayes.snp.HypothesesPrior;
import com.rtg.variant.bayes.snp.StatisticsSnp;
import com.rtg.variant.util.arithmetic.SimplePossibility;

/**
 */
public class ContaminatedSomaticCallerTest extends AbstractSomaticCallerTest<Description> {

  @Override
  protected Hypotheses<Description> getCancerHypotheses(final double same, final int ref) {
    final HypothesesPrior<Description> hyps = simpleHomoHyps(0.99, ref);
    return new HypothesesCancer(hyps, SimplePossibility.SINGLETON);
  }

  @Override
  protected AbstractSomaticCaller getSomaticCaller(final double mutation, Hypotheses<Description> hypotheses, String normalName, String cancerName, VariantParams params) {
    return new ContaminatedSomaticCaller(SomaticPriors.makeQ(mutation, 0.0, hypotheses), SomaticPriors.makeQ(mutation, 0.0, hypotheses), params);
  }

  @Override
  protected List<ModelInterface<Description>> getModel() {
    final List<ModelInterface<Description>> models = new ArrayList<>();
    for (int ref = 0; ref < 4; ref++) {
      final Hypotheses<Description> hyps = simpleHomoHyps(0.99, ref);
      final HypothesesCancer hypc = new HypothesesCancer(hyps, SimplePossibility.SINGLETON);
      models.add(new ModelCancerContamination(hypc, 0.0, new StatisticsSnp(hypc.description())));
    }
    return models;
  }
}
