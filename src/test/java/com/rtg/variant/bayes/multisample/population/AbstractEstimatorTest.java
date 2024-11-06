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

package com.rtg.variant.bayes.multisample.population;



import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import com.rtg.util.InvalidParamsException;
import com.rtg.variant.GenomePriorParams;
import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.EvidenceInterface;
import com.rtg.variant.bayes.MockEvidence;
import com.rtg.variant.bayes.MockHypotheses;
import com.rtg.variant.bayes.MockModel;
import com.rtg.variant.bayes.ModelInterface;
import com.rtg.variant.bayes.multisample.HaploidDiploidHypotheses;
import com.rtg.variant.bayes.multisample.PriorContainer;
import com.rtg.variant.bayes.multisample.population.Convergence.SimulationResult;
import com.rtg.variant.bayes.snp.DescriptionCommon;
import com.rtg.variant.bayes.snp.HypothesesNone;
import com.rtg.variant.bayes.snp.HypothesesPrior;
import com.rtg.variant.bayes.snp.HypothesesSnp;
import com.rtg.variant.bayes.snp.StatisticsSnp;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;
import com.rtg.variant.util.arithmetic.SimplePossibility;

import junit.framework.TestCase;

/**
 */
public abstract class AbstractEstimatorTest extends TestCase {

  static double[] uniform(final int size) {
    final double[] prior = new double[size];
    final double d = 1.0 / size;
    for (int i = 0; i < size; ++i) {
      prior[i] = d;
    }
    return prior;
  }

  protected abstract Estimator getEstimator();

  //Doesn't test much - just that works when all evidence the same
  public void test() {
    final double third = 1.0 / 3.0;
    final PossibilityArithmetic arith = SimplePossibility.SINGLETON;
    final DescriptionCommon descr = new DescriptionCommon("X", "Y", "Z");
    final HypothesesPrior<Description> haploid = new MockHypotheses<>(descr, arith, true, uniform(3), 0);
    final HypothesesPrior<Description> diploid = new MockHypotheses<>(descr, arith, false, uniform(6), 0);
    final double[] prob = {third, third, third};
    final List<ModelInterface<?>> models = new ArrayList<>();
    for (int i = 0; i < 100; ++i) {
      final MockModel<Description> model = new MockModel<>(diploid, new StatisticsSnp(diploid.description()), null);
      models.add(model);
      final EvidenceInterface di = new MockEvidence(descr, 0.0, prob, 1);
      for (int j = 0; j < 20; ++j) {
        model.increment(di);
      }
    }
    //final HypothesisScore[] popCalls =
    new EmAlgorithm(getEstimator(), 50).getBestScores(models, new PriorContainer<>(new HaploidDiploidHypotheses<>(HypothesesNone.SINGLETON, haploid, diploid), null));
    //System.err.println(DiploidEstimator.callToString(diploid.size(), popCalls));
  }

  //test when actual probabilities and priors agree - straightforward case -- uses human priors
  public void test1() throws InvalidParamsException {
    final PossibilityArithmetic arith = SimplePossibility.SINGLETON;
    final GenomePriorParams params = GenomePriorParams.builder().create();
    final HypothesesPrior<Description> diploid = new HypothesesSnp(arith, params, false, 0);
    final HypothesesPrior<Description> haploid = new HypothesesSnp(arith, params, true, 0);
    final List<ModelInterface<?>> models = new ArrayList<>();
    for (int i = 0; i < 100; ++i) {
      final MockModel<Description> model = new MockModel<>(diploid, new StatisticsSnp(diploid.description()), null);
      models.add(model);
    }
    final Random random = new Random(143);
    final Convergence convergence = new Convergence(haploid, diploid, getEstimator(), models, random);
    for (int i = 0; i < 5; ++i) {
      final SimulationResult iterate = convergence.simulate();
      assertEquals(0, iterate.totalIncorrect());
      //System.err.println(iterate.totalIncorrect());
    }
  }
}
