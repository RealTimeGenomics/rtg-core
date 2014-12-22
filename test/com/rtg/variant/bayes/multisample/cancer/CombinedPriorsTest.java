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

import com.rtg.util.TestUtils;
import com.rtg.util.Utils;
import com.rtg.variant.bayes.Hypotheses;
import com.rtg.variant.bayes.snp.DescriptionCommon;
import com.rtg.variant.bayes.snp.DescriptionSnp;
import com.rtg.variant.util.arithmetic.SimplePossibility;

import junit.framework.TestCase;

/**
 */
public class CombinedPriorsTest extends TestCase {

  private static class MockCombine extends CombinedPriorsSnp<DescriptionCommon> {

    private final StringBuilder mSb = new StringBuilder();

    private double mSum = 0.0;

    MockCombine(Hypotheses<DescriptionCommon> hypotheses, double mutation, double loh) {
      super(hypotheses, mutation, loh);
      integrity();
    }

    @Override
    void update(int key1, int key2, double probability) {
      mSb.append(mHypotheses.name(key1)).append(" ").append(mHypotheses.name(key2)).append(" ").append(Utils.realFormat(probability, 3)).append(LS);
      mSum += probability;
    }

    @Override
    public String toString() {
      return mSb.toString() + Utils.realFormat(mSum, 4) + LS;
    }
  }

  private static class MockHypotheses extends Hypotheses<DescriptionCommon> {
    MockHypotheses(DescriptionCommon description, boolean haploid) {
      super(description, SimplePossibility.SINGLETON, haploid);
    }
    @Override
    public int reference() {
      return 0;
    }
  }

  public void testHaploid() {
    final DescriptionCommon desc = new DescriptionCommon("A", "C");
    final Hypotheses<DescriptionCommon> hyp = new MockHypotheses(desc, true);
    final MockCombine mc = new MockCombine(hyp, 0.3, 0.1/*check this is ignored*/);
    mc.update();
    final String exp = ""
        + "A A 0.700" + LS
        + "A C 0.300" + LS
        + "C A 0.300" + LS
        + "C C 0.700" + LS
        + "2.0000" + LS
        ;
    assertEquals(exp, mc.toString());
  }


  public void testDiploidSnp() {
    final Hypotheses<DescriptionCommon> hyp = new MockHypotheses(DescriptionSnp.SINGLETON, false);
    final MockCombine mc = new MockCombine(hyp, 0.3, 0.0);
    mc.update();
    final String str = mc.toString();
    //System.err.println(str);
    TestUtils.containsAll(str,
        "10.0000",
        "A:C A:A 0.070", "A:C A:C 0.010", "A:C A:G 0.010", "A:C A:C 0.490", "A:C C:C 0.070",
        "T:T A:A 0.010", "T:T A:C 0.010", "T:T A:T 0.070", "T:T T:T 0.490"
        );
  }

  private static final String DIPLOID_EXP = ""
      + "A:A A:A 0.490" + LS
      + "A:A A:C 0.210" + LS
      + "A:A A:C 0.210" + LS
      + "A:A C:C 0.090" + LS
      + "C:C A:A 0.090" + LS
      + "C:C A:C 0.210" + LS
      + "C:C A:C 0.210" + LS
      + "C:C C:C 0.490" + LS
      + "A:C A:A 0.210" + LS
      + "A:C A:C 0.090" + LS
      + "A:C A:C 0.490" + LS
      + "A:C C:C 0.210" + LS
      ;

  public void testDiploid() {
    final DescriptionCommon desc = new DescriptionCommon("A", "C");
    final Hypotheses<DescriptionCommon> hyp = new MockHypotheses(desc, false);
    final MockCombine mc = new MockCombine(hyp, 0.3, 0.0);
    mc.update();
    final String exp = ""
        + DIPLOID_EXP
        + "3.0000" + LS
        ;
    assertEquals(exp, mc.toString());
  }

  private static final String LOH_EXP = ""
      + "A:A A:A 0.035" + LS
      + "A:A C:C 0.015" + LS
      + "A:A A:A 0.035" + LS
      + "A:A C:C 0.015" + LS
      + "C:C A:A 0.015" + LS
      + "C:C C:C 0.035" + LS
      + "C:C A:A 0.015" + LS
      + "C:C C:C 0.035" + LS
      + "A:C A:A 0.015" + LS
      + "A:C C:C 0.035" + LS
      + "A:C A:A 0.035" + LS
      + "A:C C:C 0.015" + LS
      ;

  public void testLoh() {
    final DescriptionCommon desc = new DescriptionCommon("A", "C");
    final Hypotheses<DescriptionCommon> hyp = new MockHypotheses(desc, false);
    final MockCombine mc = new MockCombine(hyp, 0.3, 0.1);
    mc.updateLoh();
    final String exp = ""
        + LOH_EXP
        + "0.3000" + LS
        ;
    assertEquals(exp, mc.toString());
  }

  public void testDiploidLoh() {
    final DescriptionCommon desc = new DescriptionCommon("A", "C");
    final Hypotheses<DescriptionCommon> hyp = new MockHypotheses(desc, false);
    final MockCombine mc = new MockCombine(hyp, 0.3, 0.1);
    mc.update();
    final String exp = ""
        //diploid * 0.9
        + "A:A A:A 0.441" + LS
        + "A:A A:C 0.189" + LS
        + "A:A A:C 0.189" + LS
        + "A:A C:C 0.081" + LS
        + "C:C A:A 0.081" + LS
        + "C:C A:C 0.189" + LS
        + "C:C A:C 0.189" + LS
        + "C:C C:C 0.441" + LS
        + "A:C A:A 0.189" + LS
        + "A:C A:C 0.081" + LS
        + "A:C A:C 0.441" + LS
        + "A:C C:C 0.189" + LS
        + LOH_EXP
        + "3.0000" + LS
        ;
    assertEquals(exp, mc.toString());
  }


}
