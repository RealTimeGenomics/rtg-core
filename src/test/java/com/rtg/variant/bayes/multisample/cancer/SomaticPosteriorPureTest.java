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

import static com.rtg.util.StringUtils.LS;

import com.rtg.variant.SomaticParams;
import com.rtg.variant.SomaticParamsBuilder;
import com.rtg.variant.VariantParams;
import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.ModelInterface;
import com.rtg.variant.bayes.snp.HypothesesPrior;

import junit.framework.TestCase;

/**
 */
public class SomaticPosteriorPureTest extends TestCase {

  protected static final String EXPECT_ALL_DIFFERENT = ""
    + "        A       C       G       T" + LS
      + "A -21.393 -29.399 -18.974 -29.399 -18.889" + LS
      + "C -24.668 -16.662 -14.243 -24.668 -14.158" + LS
      + "G -35.092 -35.092 -16.662 -35.092 -16.662" + LS
      + "T -35.092 -35.092 -24.668 -27.087 -24.583" + LS
      + "  -21.356 -16.662 -14.150 -24.574" + LS
      + "best[1,2]=-14.243" + LS
      + "equal=-15.965  notequal=-14.234" + LS
      ;

  @SuppressWarnings("unchecked")
  public <D extends Description> void testPosteriorAllDifferent() {
    final ModelInterface<Description> model = PureSomaticCallerTest.SEEN_3_C.get(0);
    final HypothesesPrior<D> hypotheses = (HypothesesPrior<D>) model.hypotheses();
    final VariantParams params = VariantParams.builder().somaticParams(getSomaticRateParams(0.001)).create();
    final AbstractSomaticCaller ccs = new PureSomaticCaller(new DefaultSomaticPriorsFactory<>(hypotheses, 0), new DefaultSomaticPriorsFactory<>(hypotheses, 0), params, 1, 1);
    ccs.integrity();

    final int length = hypotheses.size();
    final double[][] q = new double[length][length];
    new SomaticPriors<D>(hypotheses, 0.001, 0.0, DefaultSomaticPriorsFactory.defaultUniformPriors(hypotheses.description().size())) {
      @Override
      void update(int i1, int i2, double probability) {
        q[i1][i2] += probability;
      }
    }.update();

    final AbstractSomaticPosterior post = new SomaticPosteriorPure(q, model, PureSomaticCallerTest.SEEN_3_G.get(0), hypotheses, 1, 1);
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
      + "        A       C       G       T" + LS
      + "A  -0.544 -18.974 -18.974 -18.974  -0.544" + LS
      + "C -24.668 -27.087 -35.092 -35.092 -24.583" + LS
      + "G -24.668 -35.092 -27.087 -35.092 -24.583" + LS
      + "T -24.668 -35.092 -35.092 -27.087 -24.583" + LS
      + "   -0.544 -18.974 -18.974 -18.974" + LS
      + "best[0,0]=-0.544" + LS
      + "equal=-0.544  notequal=-17.872" + LS
      ;

  public void testPosteriorAllSame() {
    final HypothesesPrior<?> hypotheses = (HypothesesPrior<?>) PureSomaticCallerTest.EQUALS_REF_A.get(0).hypotheses();
    final VariantParams params = VariantParams.builder().somaticParams(getSomaticRateParams(0.001)).create();
    final AbstractSomaticCaller ccs = new PureSomaticCaller(new DefaultSomaticPriorsFactory<>(hypotheses, 0), new DefaultSomaticPriorsFactory<>(hypotheses, 0), params, 1, 1);
    ccs.integrity();

    AbstractSomaticPosterior post = null;
    for (int i = 0; i < 1; ++i) {
      final ModelInterface<Description> normal = PureSomaticCallerTest.EQUALS_REF_A.get(0);
      final ModelInterface<Description> cancer = PureSomaticCallerTest.EQUALS_REF_A.get(0);
      post = new SomaticPosteriorPure((hypotheses.haploid() ? ccs.mQHaploidFactory : ccs.mQDiploidFactory).somaticQ(0.001), normal, cancer, hypotheses, 1, 1);
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

  private static SomaticParams getSomaticRateParams(double somaticRate) {
    return new SomaticParamsBuilder().somaticRate(somaticRate).create();
  }
}
