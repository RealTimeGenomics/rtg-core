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
    final MockHypotheses<Description> hypotheses = new MockHypotheses<Description>(DescriptionSnp.SINGLETON, arith, true, priors, 0);
    final List<ModelInterface<?>> list = new ArrayList<>();
    list.add(new  MockModel<>(hypotheses, new StatisticsSnp(hypotheses.description()), null));
    final VariantParams p = VariantParams.builder()
        .machineErrorName("illumina")
        .genomePriors(GenomePriorParams.builder().create())
        .create();
    final SingletonCaller sc = new SingletonCaller(p);
    assertNull(sc.makeCall("test", 3, 4, new byte[] {0, 1, 2, 3, 4, 5}, list, new HaploidDiploidHypotheses<>(HypothesesNone.SINGLETON, null, null)));
  }

  public void testCorrect() {
    Diagnostic.setLogStream();
    final double[] priors = {0.1, 0.4, 0.35, 0.15};
    final PossibilityArithmetic arith = SimplePossibility.SINGLETON;
    final MockHypotheses<Description> hypotheses = new MockHypotheses<Description>(DescriptionSnp.SINGLETON, arith, true, priors, 2);
    final List<ModelInterface<?>> list = new ArrayList<>();
    list.add(new  Model<>(hypotheses, new StatisticsSnp(hypotheses.description()), new NoAlleleBalance()));
    final VariantParams p = VariantParams.builder()
        .machineErrorName("illumina")
        .genomePriors(GenomePriorParams.builder().create())
        .create();
    for (int i = 0; i < 10; i++) {
      list.get(0).increment(new EvidenceQ(hypotheses.description(), 3, 0, 0, 0.01, 0.1, true, false, false, false));
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
    final MockHypotheses<Description> hypotheses = new MockHypotheses<Description>(DescriptionSnp.SINGLETON, arith, true, priors, 0);
    final List<ModelInterface<?>> list = new ArrayList<>();
    list.add(new  Model<>(hypotheses, new StatisticsSnp(hypotheses.description()), new NoAlleleBalance()));
    final VariantParams p = VariantParams.builder()
        .machineErrorName("illumina")
        .genomePriors(GenomePriorParams.builder().create())
        .create();
    for (int i = 0; i < 10; i++) {
      list.get(0).increment(new EvidenceQ(hypotheses.description(), 3, 0, 0, 0.01, 0.1, true, false, false, false));
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
