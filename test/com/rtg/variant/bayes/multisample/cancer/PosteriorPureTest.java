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

import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.ModelInterface;
import com.rtg.variant.bayes.snp.HypothesesPrior;

import junit.framework.TestCase;

/**
 */
public class PosteriorPureTest extends TestCase {

  protected static final String EXPECT_ALL_DIFFERENT = ""
      + "A -21.393 -29.399 -18.974 -29.399 -18.889" + LS
      + "C -24.668 -16.662 -14.243 -24.668 -14.158" + LS
      + "G -35.092 -35.092 -16.662 -35.092 -16.662" + LS
      + "T -35.092 -35.092 -24.668 -27.087 -24.583" + LS
      + "  -21.356 -16.662 -14.150 -24.574" + LS
      + "best[1,2]=-14.243" + LS
      + "equal=-15.965  notequal=-14.234" + LS
      ;

  @SuppressWarnings("unchecked")
  public static <D extends Description> void testPosteriorAllDifferent() {
    final ModelInterface<Description> model = PureSomaticCallerTest.SEEN_3_C.get(0);
    final HypothesesPrior<D> hypotheses = (HypothesesPrior<D>) model.hypotheses();
    final AbstractSomaticCaller ccs = new PureSomaticCaller(SomaticPriors.makeQ(0.001, 0.0, hypotheses), SomaticPriors.makeQ(0.001, 0.0, hypotheses), null);
    ccs.integrity();

    final int length = hypotheses.size();
    final double[][] mQ = new double[length][length];
    new SomaticPriors<D>(hypotheses, 0.001, 0.0, SomaticPriors.defaultUniformPriors(hypotheses.description().size())) {
      @Override
      void update(int i1, int i2, double probability) {
        mQ[i1][i2] += probability;
      }
    }.update();

    final AbstractPosterior post = new PosteriorPure(mQ, model, PureSomaticCallerTest.SEEN_3_G.get(0), hypotheses);
    assertEquals(EXPECT_ALL_DIFFERENT, post.toString());
    assertEquals(1, post.bestNormal());
    assertEquals(2, post.bestCancer());
    assertEquals(1.6731, post.posteriorScore(), 1e-3);
    assertEquals(2.4018, post.normalMeasure().bestPosterior(), 2e-4);
    assertEquals(2.5032, post.cancerMeasure().bestPosterior(), 1e-4);
    assertEquals(1.731, post.ncScore(), 1e-3);
    assertFalse(post.isSameCall());
  }

  protected static final String EXPECT_ALL_SAME = ""
      + "A  -0.544 -18.974 -18.974 -18.974  -0.544" + LS
      + "C -24.668 -27.087 -35.092 -35.092 -24.583" + LS
      + "G -24.668 -35.092 -27.087 -35.092 -24.583" + LS
      + "T -24.668 -35.092 -35.092 -27.087 -24.583" + LS
      + "   -0.544 -18.974 -18.974 -18.974" + LS
      + "best[0,0]=-0.544" + LS
      + "equal=-0.544  notequal=-17.872" + LS
      ;

  public void testPosteriorAllSame() {
    HypothesesPrior<?> hypotheses = (HypothesesPrior<?>) PureSomaticCallerTest.EQUALS_REF_A.get(0).hypotheses();
    final AbstractSomaticCaller ccs = new PureSomaticCaller(SomaticPriors.makeQ(0.001, 0.0, hypotheses), SomaticPriors.makeQ(0.001, 0.0, hypotheses), null);
    ccs.integrity();

    AbstractPosterior post = null;
    for (int i = 0; i < 1; i++) {
      post = new PosteriorPure(hypotheses.haploid() ? ccs.mQHaploid : ccs.mQDiploid, PureSomaticCallerTest.EQUALS_REF_A.get(0), PureSomaticCallerTest.EQUALS_REF_A.get(0), hypotheses);
    }
    assertEquals(EXPECT_ALL_SAME, post.toString());
    assertEquals(0, post.bestNormal());
    assertEquals(0, post.bestCancer());
    assertEquals(17.32775, post.posteriorScore(), 1e-4);
    assertEquals(22.93984, post.normalMeasure().bestPosterior(), 1e-4);
    assertEquals(17.3311, post.cancerMeasure().bestPosterior(), 2e-4);
    assertEquals(-17.328, post.ncScore(), 1e-3);
    assertTrue(post.isSameCall());
  }
}
