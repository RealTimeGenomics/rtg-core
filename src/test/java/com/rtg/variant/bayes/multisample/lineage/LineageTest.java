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
package com.rtg.variant.bayes.multisample.lineage;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;

import com.rtg.variant.VariantOutputLevel;
import com.rtg.variant.VariantParams;
import com.rtg.variant.VariantSample;
import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.MockModel;
import com.rtg.variant.bayes.ModelInterface;
import com.rtg.variant.bayes.multisample.AbstractMultisampleCaller;
import com.rtg.variant.bayes.multisample.HaploidDiploidHypotheses;
import com.rtg.variant.bayes.snp.DescriptionCommon;
import com.rtg.variant.bayes.snp.HypothesesNone;
import com.rtg.variant.bayes.snp.HypothesesPrior;
import com.rtg.variant.bayes.snp.StatisticsSnp;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;
import com.rtg.variant.util.arithmetic.SimplePossibility;

import junit.framework.TestCase;

/**
 * Test class.
 */
public class LineageTest extends TestCase {

  public void testBasicLineageStructure() {
    final Lineage.LineageBuilder lb = new Lineage.LineageBuilder();
    lb.add(0, 1);
    lb.add(0, 2);
    lb.add(1, 3);
    lb.deNovoPriorDefault(0.001);
    lb.deNovoPrior(1, 0.1);
    final Lineage lineage = lb.create();
    assertEquals(0, lineage.parent(1));
    assertEquals(-1, lineage.parent(0));
    assertEquals(0, lineage.parent(2));
    assertEquals(1, lineage.parent(3));
    final HashSet<Integer> c = new HashSet<>();
    c.add(1);
    c.add(2);
    assertEquals(c, lineage.children(0));
    c.clear();
    c.add(3);
    assertEquals(c, lineage.children(1));
    assertEquals(Collections.<Integer>emptySet(), lineage.children(2));
    assertEquals(Collections.<Integer>emptySet(), lineage.children(3));
    assertTrue(lineage.isRoot(0));
    assertFalse(lineage.isRoot(2));
    assertEquals(0.1, lineage.deNovoPrior(1), 1e-10);
    assertEquals(0.001, lineage.deNovoPrior(0), 1e-10);
    assertEquals(0.001, lineage.deNovoPrior(2), 1e-10);
  }

  public void testCall() {
    final VariantParams vParams = VariantParams.builder().callLevel(VariantOutputLevel.INTERESTING).create();
    final Lineage.LineageBuilder lb = new Lineage.LineageBuilder();
    lb.params(vParams);
    lb.deNovoPriorDefault(0.0001);
    lb.add(0, 1);
    lb.add(1, 2);
    lb.add(2, 3);
    lb.add(3, 4);
    final Lineage lineage = lb.create();
    final List<ModelInterface<?>> models = new ArrayList<>();
    final DescriptionCommon descriptionCommon = new DescriptionCommon("A", "B", "C");
    final PossibilityArithmetic arithmetic = SimplePossibility.SINGLETON;
    final HypothesesPrior<Description> h = new HypothesesPrior<>(descriptionCommon, arithmetic, new double[]{0.2, 0.5, 0.3}, true, 0);
    final HypothesesPrior<Description> hypDiploid = new HypothesesPrior<>(descriptionCommon, arithmetic, new double[]{0.2, 0.5, 0.3}, false, 0);
    models.add(new MockModel<>(h.hypotheses(), new StatisticsSnp(h.description()), new double[] {0.0005, 0.0005, 0.999}));
    models.add(new MockModel<>(h.hypotheses(), new StatisticsSnp(h.description()), new double[] {0.0005, 0.0005, 0.999}));
    models.add(new MockModel<>(h.hypotheses(), new StatisticsSnp(h.description()), new double[] {0.999, 0.0005, 0.0005}));
    models.add(new MockModel<>(h.hypotheses(), new StatisticsSnp(h.description()), new double[] {0.999, 0.0009, 0.0001}));
    models.add(new MockModel<>(h.hypotheses(), new StatisticsSnp(h.description()), new double[] {0.999, 0.0009, 0.0001}));
    final HaploidDiploidHypotheses<HypothesesPrior<Description>> haploidDiploid = new HaploidDiploidHypotheses<>(HypothesesNone.SINGLETON, h, hypDiploid);
    final AbstractMultisampleCaller.ComparisonResult comparisonResult = lineage.makeSamples(models, haploidDiploid);
    assertNotNull(comparisonResult);
    final String[] strings = new String[]{"C", "C", "A", "A", "A"};
    final VariantSample.DeNovoStatus[] status = {VariantSample.DeNovoStatus.UNSPECIFIED, VariantSample.DeNovoStatus.NOT_DE_NOVO, VariantSample.DeNovoStatus.IS_DE_NOVO, VariantSample.DeNovoStatus.NOT_DE_NOVO, VariantSample.DeNovoStatus.NOT_DE_NOVO};
    assertTrue(comparisonResult.isInteresting());
    for (int i = 0; i < strings.length; ++i) {
      final String expected = strings[i];
      final VariantSample variantSample = comparisonResult.getSamples()[i];
      assertEquals("Sample: " + i, expected, variantSample.getName());
      assertEquals("Sample: " + i, status[i], variantSample.isDeNovo());
    }
  }

  public void testUninteresting() {
    final VariantParams vParams = VariantParams.builder().callLevel(VariantOutputLevel.INTERESTING).create();
    final Lineage.LineageBuilder lb = new Lineage.LineageBuilder();
    lb.params(vParams);
    lb.deNovoPriorDefault(0.0001);
    lb.add(0, 1);
    lb.add(1, 2);
    lb.add(2, 3);
    lb.add(3, 4);
    final Lineage lineage = lb.create();
    final List<ModelInterface<?>> models = new ArrayList<>();
    final DescriptionCommon descriptionCommon = new DescriptionCommon("A", "B", "C");
    final PossibilityArithmetic arithmetic = SimplePossibility.SINGLETON;
    final HypothesesPrior<Description> h = new HypothesesPrior<>(descriptionCommon, arithmetic, new double[]{0.2, 0.5, 0.3}, true, 0);
    final HypothesesPrior<Description> hypDiploid = new HypothesesPrior<>(descriptionCommon, arithmetic, new double[]{0.2, 0.5, 0.3}, false, 0);
    models.add(new MockModel<>(h.hypotheses(), new StatisticsSnp(h.description()), new double[] {0.999, 0.0005, 0.0005}));
    models.add(new MockModel<>(h.hypotheses(), new StatisticsSnp(h.description()), new double[] {0.999, 0.0005, 0.0005}));
    models.add(new MockModel<>(h.hypotheses(), new StatisticsSnp(h.description()), new double[] {0.999, 0.0005, 0.0005}));
    models.add(new MockModel<>(h.hypotheses(), new StatisticsSnp(h.description()), new double[] {0.999, 0.0009, 0.0001}));
    models.add(new MockModel<>(h.hypotheses(), new StatisticsSnp(h.description()), new double[] {0.999, 0.0009, 0.0001}));
    final HaploidDiploidHypotheses<HypothesesPrior<Description>> haploidDiploid = new HaploidDiploidHypotheses<>(HypothesesNone.SINGLETON, h, hypDiploid);
    final AbstractMultisampleCaller.ComparisonResult comparisonResult = lineage.makeSamples(models, haploidDiploid);
    assertNull(comparisonResult);
  }
}
