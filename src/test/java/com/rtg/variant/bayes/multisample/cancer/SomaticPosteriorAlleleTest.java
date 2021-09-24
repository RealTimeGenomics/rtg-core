/*
 * Copyright (c) 2017. Real Time Genomics Limited.
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
