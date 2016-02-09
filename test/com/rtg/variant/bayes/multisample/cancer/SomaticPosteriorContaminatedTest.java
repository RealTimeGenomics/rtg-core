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

import static com.rtg.util.StringUtils.LS;

import com.rtg.variant.SomaticParamsBuilder;
import com.rtg.variant.VariantParams;
import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.Evidence;
import com.rtg.variant.bayes.Hypotheses;
import com.rtg.variant.bayes.Model;
import com.rtg.variant.bayes.ModelInterface;
import com.rtg.variant.bayes.NoAlleleBalance;
import com.rtg.variant.bayes.snp.EvidenceQ;
import com.rtg.variant.bayes.snp.HypothesesPrior;
import com.rtg.variant.bayes.snp.StatisticsSnp;
import com.rtg.variant.dna.DNARange;
import com.rtg.variant.util.arithmetic.SimplePossibility;

import junit.framework.TestCase;

/**
 */
public class SomaticPosteriorContaminatedTest extends TestCase {

  public void testPosteriorAllDifferent00() {
    final double contamination = 0.0;
    final AbstractSomaticPosterior post = getContaminatedPosterior(contamination);
    assertEquals(SomaticPosteriorPureTest.EXPECT_ALL_DIFFERENT, post.toString());
    assertEquals(1, post.bestNormal());
    assertEquals(2, post.bestCancer());
    assertEquals(1.6731, post.posteriorScore(), 1e-4);
    assertEquals(2.4018, post.normalMeasure().bestPosterior(), 2e-4);
    assertEquals(2.5032, post.cancerMeasure().bestPosterior(), 1e-4);
    assertEquals(1.731, post.ncScore(), 1e-3);
    assertFalse(post.isSameCall());
  }

  private AbstractSomaticPosterior getContaminatedPosterior(double contamination) {
    final HypothesesPrior<?> hypotheses = (HypothesesPrior<?>) PureSomaticCallerTest.SEEN_3_C.get(0).hypotheses();
    final VariantParams params = VariantParams.builder().somaticParams(new SomaticParamsBuilder().somaticRate(0.001).create()).create();
    final ContaminatedSomaticCaller cc = new ContaminatedSomaticCaller(new SomaticPriorsFactory<>(hypotheses, 0), new SomaticPriorsFactory<>(hypotheses, 0), params, 1, 1);
    cc.integrity();
    // construct a contaminated cancer bayesian.
    final int numReads = 3;
    final int refNt = DNARange.A;
    final int refCode = refNt - 1;
    final Hypotheses<Description> simpleHomoHyps = AbstractSomaticCallerTest.simpleHomoHyps(0.99, refCode);
    final HypothesesCancer<Hypotheses<Description>> hypc = new HypothesesCancer<>(simpleHomoHyps, SimplePossibility.SINGLETON);
    final ModelCancerContamination<Hypotheses<Description>> cancer = new ModelCancerContamination<>(hypc, contamination, new StatisticsSnp(hypc.description()), new NoAlleleBalance());
    final ModelInterface<Description> normal = new Model<>(simpleHomoHyps, new StatisticsSnp(simpleHomoHyps.description()), new NoAlleleBalance());
    // run through several identical reads
    final Evidence evc = new EvidenceQ(simpleHomoHyps.description(), 1, 0, 0, 0.05, 0.05, true, false, false, false);
    final Evidence evg = new EvidenceQ(simpleHomoHyps.description(), 2, 0, 0, 0.05, 0.05, true, false, false, false);

    for (int i = 0; i < numReads; i++) {
      cancer.increment(evg);
      normal.increment(evc);
    }
    cancer.freeze();
    normal.freeze();
    return new SomaticPosteriorContaminated((normal.haploid() ? cc.mQHaploidFactory : cc.mQDiploidFactory).somaticQ(0.001), normal, cancer, hypotheses, 1, 1);
  }

  private static final String EXPECT_ALL_DIFFERENT_20 = ""
      + "A -21.393 -29.399 -19.620 -29.399 -19.463" + LS
      + "C -24.668 -16.662 -14.889 -24.668 -14.732" + LS
      + "G -29.146 -29.146 -16.662 -29.146 -16.662" + LS
      + "T -35.092 -35.092 -25.314 -27.087 -25.157" + LS
      + "  -21.356 -16.662 -14.725 -24.564" + LS
      + "best[1,2]=-14.889" + LS
      + "equal=-15.965  notequal=-14.881" + LS
      ;
  public void testPosteriorAllDifferent20() {
    final double contamination = 0.20;
    final AbstractSomaticPosterior post = getContaminatedPosterior(contamination);
    assertEquals(EXPECT_ALL_DIFFERENT_20, post.toString());
    assertEquals(1, post.bestNormal());
    assertEquals(2, post.bestCancer());
    assertEquals(1.87076, post.normalMeasure().bestPosterior(), 1e-4);
    assertEquals(1.92789, post.cancerMeasure().bestPosterior(), 1e-4);
    assertEquals(1.04949, post.posteriorScore(), 1e-4);
    assertEquals(1.084, post.ncScore(), 1e-3);
    assertFalse(post.isSameCall());
  }

  public void testPosteriorAllSame() {
    final ModelInterface<Description> normal = PureSomaticCallerTest.EQUALS_REF_A.get(0);
    final HypothesesPrior<?> hypotheses = (HypothesesPrior<?>) normal.hypotheses();
    final VariantParams params = VariantParams.builder().somaticParams(new SomaticParamsBuilder().somaticRate(0.001).create()).create();
    final ContaminatedSomaticCaller cc = new ContaminatedSomaticCaller(new SomaticPriorsFactory<>(hypotheses, 0), new SomaticPriorsFactory<>(hypotheses, 0), params, 1, 1);
    cc.integrity();
    final int numReads = 3;
    final int refNt = DNARange.A;
    final int refCode = refNt - 1;
    final Hypotheses<Description> simpleHomoHyps = AbstractSomaticCallerTest.simpleHomoHyps(0.99, refCode);
    final HypothesesCancer<Hypotheses<Description>> hypc = new HypothesesCancer<>(simpleHomoHyps, SimplePossibility.SINGLETON);
    final ModelCancerContamination<Hypotheses<Description>> cancer = new ModelCancerContamination<>(hypc, 0.0, new StatisticsSnp(hypc.description()), new NoAlleleBalance());
    final Evidence eva = new EvidenceQ(simpleHomoHyps.description(), 0, 0, 0, 0.05, 0.05, true, false, false, false);
    // run through several identical reads
    for (int i = 0; i < numReads; i++) {
      cancer.increment(eva);
    }
    cancer.freeze();
    final AbstractSomaticPosterior post = new SomaticPosteriorContaminated((hypotheses.haploid() ? cc.mQHaploidFactory : cc.mQDiploidFactory).somaticQ(0.001), normal, cancer, hypotheses, 1, 1);
    assertEquals(SomaticPosteriorPureTest.EXPECT_ALL_SAME, post.toString());
    assertEquals(0, post.bestNormal());
    assertEquals(0, post.bestCancer());
    assertEquals(17.32775, post.posteriorScore(), 1e-4);
    assertEquals(22.93984, post.normalMeasure().bestPosterior(), 1e-4);
    assertEquals(17.3311, post.cancerMeasure().bestPosterior(), 1e-4);
    assertEquals(-17.328, post.ncScore(), 1e-3);
    assertTrue(post.isSameCall());
  }

}
