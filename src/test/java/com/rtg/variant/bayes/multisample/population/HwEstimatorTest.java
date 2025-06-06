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

import com.rtg.variant.GenomePriorParams;
import com.rtg.variant.VariantSample;
import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.MockModel;
import com.rtg.variant.bayes.ModelInterface;
import com.rtg.variant.bayes.multisample.HaploidDiploidHypotheses;
import com.rtg.variant.bayes.multisample.HypothesisScore;
import com.rtg.variant.bayes.multisample.PriorContainer;
import com.rtg.variant.bayes.snp.DescriptionCommon;
import com.rtg.variant.bayes.snp.EvidenceQ;
import com.rtg.variant.bayes.snp.HypothesesCommon;
import com.rtg.variant.bayes.snp.HypothesesNone;
import com.rtg.variant.bayes.snp.HypothesesPrior;
import com.rtg.variant.bayes.snp.HypothesesSnp;
import com.rtg.variant.bayes.snp.StatisticsSnp;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;
import com.rtg.variant.util.arithmetic.SimplePossibility;

/**
 */
public class HwEstimatorTest extends AbstractEstimatorTest {

  private HaploidDiploidHypotheses<HypothesesPrior<Description>> mHDH;
  private List<ModelInterface<?>> mModels;

  @Override
  protected Estimator getEstimator() {
    return new HwEstimator();
  }

  @Override
  public void test() {
    assertEquals("HwEstimator", getEstimator().toString());
  }

  @Override
  public void setUp() throws Exception {
    super.setUp();

    final PossibilityArithmetic arith = SimplePossibility.SINGLETON;
    final GenomePriorParams params = GenomePriorParams.builder().create();

    final Description d = new DescriptionCommon("A", "T");
    final HypothesesPrior<Description> diploid = new HypothesesCommon<Description>(d, arith, false, 0) {
      @Override
      public double p(int hyp) {
        return 0.3333;
      }
    };
    final HypothesesSnp haploid = new HypothesesSnp(arith, params, true, 0);

    final DescriptionCounts dc = new DescriptionCounts(haploid.size(), 0);
    dc.increment(0, 2);
    dc.increment(1, 4);
    dc.increment(2, 1);

    mHDH = new HaploidDiploidHypotheses<>(HypothesesNone.SINGLETON, haploid, diploid, true, dc);
    mModels = new ArrayList<>();

    final MockModel<Description> model = new MockModel<>(diploid, new StatisticsSnp(diploid.description()), null);
    for (int i = 0; i < 6; ++i) {
      model.increment(new EvidenceQ(d, 1, 10, 10, 0.0, 0.1, true, true, true, true, false));
    }
    for (int i = 0; i < 7; ++i) {
      model.increment(new EvidenceQ(d, 0, 10, 10, 0.0, 0.1, true, true, true, true, false));
    }
    mModels.add(model);
    final MockModel<Description> model2 = new MockModel<>(haploid, new StatisticsSnp(haploid.description()), null);
    mModels.add(model2);

    final MockModel<Description> denovoModel = new MockModel<Description>(haploid, new StatisticsSnp(haploid.description()), null) {
      @Override
      public HypothesisScore best(HypothesesPrior<?> hypotheses) {
        final HypothesisScore hs = super.best(hypotheses);
        hs.setDenovo(VariantSample.DeNovoStatus.IS_DE_NOVO);
        return hs;
      }
    };
    mModels.add(denovoModel);

    final MockModel<Description> nullModel = new MockModel<Description>(haploid, new StatisticsSnp(haploid.description()), null) {
      @Override
      public HypothesisScore best(HypothesesPrior<?> hypotheses) {
        return null;
      }
    };
    mModels.add(nullModel);
    for (ModelInterface<?> modelInterface : mModels) {
      modelInterface.freeze();
    }
  }


  @Override
  public void tearDown() throws Exception {
    super.tearDown();
    mHDH = null;
    mModels = null;
  }

  public void testEstimate() {
    final HwEstimator e = new HwEstimator() {

      @Override
      public <D extends Description, T extends HypothesesPrior<D>> HaploidDiploidHypotheses<HypothesesPrior<D>> computeNewPriors(HaploidDiploidHypotheses<T> hypotheses, int[] haploidCounts, int haploidTotal) {

        assertNotNull(haploidCounts);
        assertEquals(10, haploidTotal);

        assertEquals(4, haploidCounts.length);

        assertEquals(4, haploidCounts[0]); //5
        assertEquals(5, haploidCounts[1]); //4
        assertEquals(1, haploidCounts[2]);
        assertEquals(0, haploidCounts[3]);

        assertEquals(mHDH, hypotheses);

        return null;
      }
    };

    final PriorContainer<HypothesesPrior<Description>> pc = new PriorContainer<>(mHDH, null);

    final EmResult<?> res = e.estimate(mModels, pc);

    assertNotNull(res);
    assertNotNull(res.getPriorContainer());
    assertNull(res.getPriorContainer().getHypotheses());
    assertNotNull(res.calls());
    assertNull(res.getPriorContainer().getBs());
  }

  public void testComputePriors() {
    final HwEstimator hw = new HwEstimator();

    final PriorContainer<HypothesesPrior<Description>> pc = new PriorContainer<>(mHDH, null);

    final EmResult<?> res = hw.estimate(mModels, pc);

    assertNotNull(res);
    assertNotNull(res.getPriorContainer());
    assertNotNull(res.calls());

    final HaploidDiploidHypotheses<? extends HypothesesPrior<?>> hdh = res.getPriorContainer().getHypotheses();

    assertNotNull(hdh);
    assertNotSame(mHDH, hdh);
    assertFalse(hdh.isDefault());

    final double delta = 1e-3;

    final HypothesesPrior<? extends Description> hap = hdh.haploid();
    assertNotNull(hap);
    assertEquals(4, hap.size());
    assertEquals(0.454, hap.p(0), delta);
    assertEquals(0.454, hap.p(1), delta);
    assertEquals(0.090, hap.p(2), delta);
    assertEquals(5.11E-6, hap.p(3), delta);

    final HypothesesPrior<?> diploid = hdh.diploid();
    assertNotNull(diploid);

    assertEquals(3, diploid.size());
    assertEquals(0.201, diploid.p(0), delta);
    assertEquals(0.209, diploid.p(1), delta);
    assertEquals(0.415, diploid.p(2), delta);
  }

  public void testOverflow() {
    final HwEstimator hw = new HwEstimator();
    final DescriptionCommon d = new DescriptionCommon("A", "C", "G", "T");
    final HypothesesCommon<DescriptionCommon> hap = new HypothesesCommon<DescriptionCommon>(d, SimplePossibility.SINGLETON, true, 0) {
      @Override
      public double p(int hyp) {
        return 0.25;
      }
    };
    final HypothesesCommon<DescriptionCommon> dip = new HypothesesCommon<DescriptionCommon>(d, SimplePossibility.SINGLETON, false, 0) {
      @Override
      public double p(int hyp) {
        return 0.25;
      }
    };
    final HaploidDiploidHypotheses<HypothesesCommon<DescriptionCommon>> hdh = new HaploidDiploidHypotheses<>(null, hap, dip);
    hw.computeNewPriors(hdh, new int[] {0, 1, 111581, 0}, 111583);
  }
}
