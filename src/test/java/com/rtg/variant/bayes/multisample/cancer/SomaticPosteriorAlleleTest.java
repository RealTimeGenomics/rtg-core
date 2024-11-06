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

package com.rtg.variant.bayes.multisample.cancer;

import com.rtg.variant.SomaticParamsBuilder;
import com.rtg.variant.VariantParams;
import com.rtg.variant.bayes.HypothesesPowerSet;
import com.rtg.variant.bayes.MockEvidence;
import com.rtg.variant.bayes.Model;
import com.rtg.variant.bayes.NoAlleleBalance;
import com.rtg.variant.bayes.StatisticsDouble;
import com.rtg.variant.bayes.snp.DescriptionCommon;
import com.rtg.variant.bayes.snp.HypothesesCommon;
import com.rtg.variant.bayes.snp.HypothesesPrior;
import com.rtg.variant.util.arithmetic.LogPossibility;

import junit.framework.TestCase;

/**
 * Test the corresponding class.
 */
public class SomaticPosteriorAlleleTest extends TestCase {

  private static final double[][] P = {
    {0.97, 0.01, 0.01, 0.01},
    {0.01, 0.97, 0.01, 0.01},
    {0.01, 0.01, 0.97, 0.01},
    {0.01, 0.01, 0.01, 0.97},
  };

  private void increment(final Model<?> model, final int base) {
    model.increment(new MockEvidence(model.description(), 0.05, P[base], base));
  }

  public void testPosteriorHaploidNormal() {
    final double mu = 0.01;
    final DescriptionCommon desc = new DescriptionCommon("A", "C", "G", "T");
    final HypothesesCommon<DescriptionCommon> normalHyp = new HypothesesCommon<>(desc, LogPossibility.SINGLETON, true, 0);
    final HypothesesPowerSet<DescriptionCommon> cancerHyp = new HypothesesPowerSet<>(desc, LogPossibility.SINGLETON, 0);
    final VariantParams params = VariantParams.builder().somaticParams(new SomaticParamsBuilder().somaticRate(mu).create()).create();
    final AlleleSomaticCaller ccs = new AlleleSomaticCaller(new AlleleSomaticPriorsFactory<>(normalHyp), null, params, 1, 1);
    ccs.integrity();
    final Model<DescriptionCommon> normalModel = new Model<>(normalHyp, new StatisticsDouble(desc), new NoAlleleBalance());
    final ModelCancerAllele<DescriptionCommon> cancerModel = new ModelCancerAllele<>(cancerHyp, new StatisticsDouble(desc));
    final double[][] q = SomaticPriorsAllele.makeQ(mu, normalHyp, cancerHyp);
    for (int k = 0; k < 10; ++k) {
      increment(normalModel, 1); // C
      increment(cancerModel, 2); // G
    }
    normalModel.freeze();
    cancerModel.freeze();
    final AbstractSomaticPosterior post = new SomaticPosteriorAllele(q, normalModel, cancerModel, (HypothesesPrior<?>) normalModel.hypotheses(), 1, 1);
    post.postConstruction();
    assertEquals(1, post.bestNormal());
    //System.out.println(post.toString());
    assertEquals((1 << 2) - 1, post.bestCancer());
    assertEquals(-23.423, post.mPosterior[post.bestNormal()][post.bestCancer()], 1e-3);
  }

  public void testPosteriorDiploidNormal() {
    final double mu = 0.01;
    final DescriptionCommon desc = new DescriptionCommon("A", "C", "G", "T");
    final HypothesesCommon<DescriptionCommon> normalHyp = new HypothesesCommon<>(desc, LogPossibility.SINGLETON, false, 0);
    final HypothesesPowerSet<DescriptionCommon> cancerHyp = new HypothesesPowerSet<>(desc, LogPossibility.SINGLETON, 0);
    final VariantParams params = VariantParams.builder().somaticParams(new SomaticParamsBuilder().somaticRate(mu).create()).create();
    final AlleleSomaticCaller ccs = new AlleleSomaticCaller(null, new AlleleSomaticPriorsFactory<>(normalHyp), params, 1, 1);
    ccs.integrity();
    final Model<DescriptionCommon> normalModel = new Model<>(normalHyp, new StatisticsDouble(desc), new NoAlleleBalance());
    final ModelCancerAllele<DescriptionCommon> cancerModel = new ModelCancerAllele<>(cancerHyp, new StatisticsDouble(desc));
    final double[][] q = SomaticPriorsAllele.makeQ(mu, normalHyp, cancerHyp);
    for (int k = 0; k < 10; ++k) {
      increment(normalModel, 1); // C
      increment(normalModel, 3); // T
      increment(cancerModel, 2); // G
      increment(cancerModel, 3); // T
    }
    normalModel.freeze();
    cancerModel.freeze();
    final AbstractSomaticPosterior post = new SomaticPosteriorAllele(q, normalModel, cancerModel, (HypothesesPrior<?>) normalModel.hypotheses(), 1, 1);
    post.postConstruction();
    assertEquals(8, post.bestNormal()); // C:T
    assertEquals(((1 << 2) | (1 << 3)) - 1, post.bestCancer()); // G:T
    assertEquals(-63.933, post.mPosterior[post.bestNormal()][post.bestCancer()], 1e-3);
  }

  // These ones are for generating posterior tables of how the calling looks
//  public void testGraphAlleleSomaticCaller() {
//    final double mu = 0.000001;
//    final DescriptionCommon desc = new DescriptionCommon("A", "C", "G", "T");
//    final HypothesesCommon<DescriptionCommon> normalHyp = new HypothesesCommon<>(desc, LogPossibility.SINGLETON, false, 0);
//    final HypothesesPowerSet<DescriptionCommon> cancerHyp = new HypothesesPowerSet<>(desc, LogPossibility.SINGLETON, 0);
//    final double[][] q = SomaticPriorsAllele.makeQ(mu, normalHyp, cancerHyp);
//    for (int t = 0; t < 500; ++t) {
//      final Model<DescriptionCommon> normalModel = new Model<>(normalHyp, new StatisticsDouble(desc), new NoAlleleBalance());
//      final ModelCancerAllele<DescriptionCommon> cancerModel = new ModelCancerAllele<>(cancerHyp, new StatisticsDouble(desc));
//      for (int k = 0; k < 500; ++k) {
//        increment(normalModel, 1); // C
//        increment(cancerModel, k >= t ? 1 : 3); // C or T
//      }
//      normalModel.freeze();
//      cancerModel.freeze();
//      final AbstractSomaticPosterior post = new SomaticPosteriorAllele(q, normalModel, cancerModel, (HypothesesPrior<?>) normalModel.hypotheses(), 1, 1);
//      post.postConstruction();
//      // C, CT, T
//      System.out.println(t + " " + normalHyp.name(post.bestNormal()) + " " + cancerHyp.name(post.bestCancer()) + " " + post.mPosterior[1][1] + " " + post.mPosterior[1][9] + " " + post.mPosterior[1][7] + " " + post.mPosterior[post.bestNormal()][post.bestCancer()]);
//    }
//  }
//
//  public void testGraphContamCaller() {
//    final double mu = 0.000001;
//    final double contam = 0.75;
//    final DescriptionCommon desc = new DescriptionCommon("A", "C", "G", "T");
//    final HypothesesCommon<DescriptionCommon> normalHyp = new HypothesesCommon<>(desc, LogPossibility.SINGLETON, false, 0);
//    final HypothesesCancer<?> cancerHyp = new HypothesesCancer<>(normalHyp, LogPossibility.SINGLETON);
//    final double[][] q = new DefaultSomaticPriorsFactory<>(normalHyp, 0).somaticQ(mu);
//    for (int t = 0; t < 500; ++t) {
//      final Model<DescriptionCommon> normalModel = new Model<>(normalHyp, new StatisticsDouble(desc), new NoAlleleBalance());
//      final ModelCancerContamination<?> cancerModel = new ModelCancerContamination<>(cancerHyp, contam, new StatisticsDouble(desc), new NoAlleleBalance());
//      for (int k = 0; k < 500; ++k) {
//        increment(normalModel, 1); // C
//        increment(cancerModel, k >= t ? 1 : 3); // C or T
//      }
//      normalModel.freeze();
//      cancerModel.freeze();
//      final AbstractSomaticPosterior post = new SomaticPosteriorContaminated(q, normalModel, cancerModel, (HypothesesPrior<?>) normalModel.hypotheses(), 1, 1, contam, false);
//      post.postConstruction();
//      // C, CT, T
//      System.out.println(t + " " + normalHyp.name(post.bestNormal()) + " " + cancerHyp.name(post.bestCancer()) + " " + post.mPosterior[1][1] + " " + post.mPosterior[1][normalHyp.code().code(1, 3)] + " " + post.mPosterior[1][3] + " " + post.mPosterior[post.bestNormal()][post.bestCancer()]);
//    }
//  }
}
