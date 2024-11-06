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

package com.rtg.variant.bayes.multisample.singleton;


import java.util.ArrayList;
import java.util.List;

import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.variant.GenomePriorParams;
import com.rtg.variant.Variant;
import com.rtg.variant.VariantParams;
import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.MockHypotheses;
import com.rtg.variant.bayes.MockModel;
import com.rtg.variant.bayes.Model;
import com.rtg.variant.bayes.ModelInterface;
import com.rtg.variant.bayes.NoAlleleBalance;
import com.rtg.variant.bayes.multisample.HaploidDiploidHypotheses;
import com.rtg.variant.bayes.snp.DescriptionSnp;
import com.rtg.variant.bayes.snp.EvidenceQ;
import com.rtg.variant.bayes.snp.HypothesesNone;
import com.rtg.variant.bayes.snp.StatisticsSnp;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;
import com.rtg.variant.util.arithmetic.SimplePossibility;

import junit.framework.TestCase;

/**
 */
public class SingletonCallerTest extends TestCase {


  private static final double[] PRIORS = new double[]{0.1, 0.4, 0.35, 0.15};

  public void test() {
    final double[] priors = PRIORS;
    final PossibilityArithmetic arith = SimplePossibility.SINGLETON;
    final MockHypotheses<Description> hypotheses = new MockHypotheses<>(DescriptionSnp.SINGLETON, arith, true, priors, 0);
    final List<ModelInterface<?>> list = new ArrayList<>();
    list.add(new  MockModel<>(hypotheses, new StatisticsSnp(hypotheses.description()), null));
    final VariantParams p = VariantParams.builder()
        .machineErrorName("illumina")
        .genomePriors(GenomePriorParams.builder().create())
        .create();
    final SingletonCaller sc = new SingletonCaller(p);
    for (ModelInterface<?> modelInterface : list) {
      modelInterface.freeze();
    }
    final Variant test = sc.makeCall("test", 3, 4, new byte[]{0, 1, 2, 3, 1, 1}, list, new HaploidDiploidHypotheses<>(HypothesesNone.SINGLETON, hypotheses, null));
    assertEquals("C", test.getSample(0).getName());
  }

  public void testCorrect() {
    Diagnostic.setLogStream();
    final double[] priors = {0.1, 0.4, 0.35, 0.15};
    final PossibilityArithmetic arith = SimplePossibility.SINGLETON;
    final MockHypotheses<Description> hypotheses = new MockHypotheses<>(DescriptionSnp.SINGLETON, arith, true, priors, 2);
    final List<ModelInterface<?>> list = new ArrayList<>();
    list.add(new  Model<>(hypotheses, new StatisticsSnp(hypotheses.description()), new NoAlleleBalance()));
    final VariantParams p = VariantParams.builder()
        .machineErrorName("illumina")
        .genomePriors(GenomePriorParams.builder().create())
        .create();
    for (int i = 0; i < 10; ++i) {
      list.get(0).increment(new EvidenceQ(hypotheses.description(), 3, 0, 0, 0.01, 0.1, true, false, true, false, false));
    }
    for (ModelInterface<?> modelInterface : list) {
      modelInterface.freeze();
    }
    final SingletonCaller sc = new SingletonCaller(p);
    final Variant v = sc.makeCall("test", 3, 4, new byte[] {0, 1, 2, 3, 4, 5}, list, new HaploidDiploidHypotheses<>(HypothesesNone.SINGLETON, hypotheses, null));
    assertNotNull(v);
    assertEquals("test", v.getLocus().getSequenceName());
    //System.err.println(mvc);
    assertEquals("G", v.getLocus().getRefNts());
    assertEquals('C', v.getLocus().getPreviousRefNt());
    assertEquals(3, v.getLocus().getStart());
  }

  public void testCorrectLeftEndTemplate() {
    final double[] priors = {0.1, 0.4, 0.35, 0.15};
    final PossibilityArithmetic arith = SimplePossibility.SINGLETON;
    final MockHypotheses<Description> hypotheses = new MockHypotheses<>(DescriptionSnp.SINGLETON, arith, true, priors, 0);
    final List<ModelInterface<?>> list = new ArrayList<>();
    list.add(new  Model<>(hypotheses, new StatisticsSnp(hypotheses.description()), new NoAlleleBalance()));
    final VariantParams p = VariantParams.builder()
        .machineErrorName("illumina")
        .genomePriors(GenomePriorParams.builder().create())
        .create();
    for (int i = 0; i < 10; ++i) {
      list.get(0).increment(new EvidenceQ(hypotheses.description(), 3, 0, 0, 0.01, 0.1, true, false, true, false, false));
    }
    for (ModelInterface<?> modelInterface : list) {
      modelInterface.freeze();
    }
    final SingletonCaller sc = new SingletonCaller(p);
    final Variant v = sc.makeCall("test", 0, 1, new byte[] {1, 1, 2, 3, 4, 5}, list, new HaploidDiploidHypotheses<>(HypothesesNone.SINGLETON, hypotheses, null));
    assertNotNull(v);
    assertEquals("test", v.getLocus().getSequenceName());
    assertEquals("A", v.getLocus().getRefNts());
    assertEquals('N', v.getLocus().getPreviousRefNt());
    assertEquals(0, v.getLocus().getStart());
  }
}
