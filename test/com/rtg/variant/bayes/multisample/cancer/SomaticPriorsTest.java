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

import java.io.IOException;

import com.rtg.util.InvalidParamsException;
import com.rtg.util.TestUtils;
import com.rtg.util.Utils;
import com.rtg.util.integrity.Exam;
import com.rtg.variant.bayes.Hypotheses;
import com.rtg.variant.bayes.MockHypotheses;
import com.rtg.variant.bayes.snp.DescriptionCommon;
import com.rtg.variant.bayes.snp.DescriptionSnp;
import com.rtg.variant.util.arithmetic.SimplePossibility;

import junit.framework.TestCase;

/**
 */
public class SomaticPriorsTest extends TestCase {

  private static class MockSomaticPriors extends SomaticPriors<DescriptionCommon> {

    private final StringBuilder mSb = new StringBuilder();

    private double mSum = 0.0;

    MockSomaticPriors(Hypotheses<DescriptionCommon> hypotheses, double mutation, double loh) {
      super(hypotheses, mutation, loh, DefaultSomaticPriorsFactory.defaultUniformPriors(hypotheses.description().size()));
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

  private static class MyMockHypotheses extends Hypotheses<DescriptionCommon> {
    MyMockHypotheses(DescriptionCommon description, boolean haploid) {
      super(description, SimplePossibility.SINGLETON, haploid);
    }
    @Override
    public int reference() {
      return 0;
    }
  }

  public void testHaploidMock() {
    final DescriptionCommon desc = new DescriptionCommon("A", "C");
    final Hypotheses<DescriptionCommon> hyp = new MyMockHypotheses(desc, true);
    final MockSomaticPriors mc = new MockSomaticPriors(hyp, 0.3, 0.1/*check this is ignored*/);
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
    final Hypotheses<DescriptionCommon> hyp = new MyMockHypotheses(DescriptionSnp.SINGLETON, false);
    final MockSomaticPriors mc = new MockSomaticPriors(hyp, 0.3, 0.0);
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

  public void testDiploidMock() {
    final DescriptionCommon desc = new DescriptionCommon("A", "C");
    final Hypotheses<DescriptionCommon> hyp = new MyMockHypotheses(desc, false);
    final MockSomaticPriors mc = new MockSomaticPriors(hyp, 0.3, 0.0);
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
    final Hypotheses<DescriptionCommon> hyp = new MyMockHypotheses(desc, false);
    final MockSomaticPriors mc = new MockSomaticPriors(hyp, 0.3, 0.1);
    mc.updateLoh();
    final String exp = ""
        + LOH_EXP
        + "0.3000" + LS
        ;
    assertEquals(exp, mc.toString());
  }

  public void testDiploidLoh() {
    final DescriptionCommon desc = new DescriptionCommon("A", "C");
    final Hypotheses<DescriptionCommon> hyp = new MyMockHypotheses(desc, false);
    final MockSomaticPriors mc = new MockSomaticPriors(hyp, 0.3, 0.1);
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

  public void testMutationNormalize() {
    final double third = 1.0 / 3.0;
    final double sixth = third / 2.0;
    checkMN(0.0, 1, new double[] {third, third, third}, new double[] {0.0, 1.0, 0.0});
    checkMN(0.5, 2, new double[] {third, third, third}, new double[] {sixth, sixth, 2.0 * third});
  }

  private void checkMN(final double mutation, final int ref, final double[] prior, final double[] exp) {
    final int len = exp.length;
    assertEquals(len, prior.length);
    final double[] norm = SomaticPriors.mutationNormalize(mutation, ref, prior);
    assertEquals(len, norm.length);
    double sum = 0.0;
    for (int i = 0; i < len; ++i) {
      final double no = norm[i];
      assertEquals(exp[i], no, 1e-7);
      assertTrue(Double.isFinite(no));
      sum += no;
    }
    assertEquals(1.0, sum, 1e-7);
  }

  public void testHaploid() throws InvalidParamsException, IOException {
    final DescriptionCommon desc = new DescriptionCommon("", "A", "AA");
    final Hypotheses<DescriptionCommon> hyp = new MockHypotheses<>(desc, SimplePossibility.SINGLETON, true, new double[] {0.0, 0.0, 0.0}, 0);
    final double[][] initialPriors = {
      {0.5, 0.3, 0.2},
      {0.25, 0.5, 0.25},
      {0.5 / 3.0, 1.0 / 3.0, 0.5},
    };
    final double[][] q = SomaticPriors.makeQ(0.5, 0.0, hyp, initialPriors);
    final int n = hyp.size();
    assertEquals(n, q.length);
    //System.err.println(IntegralAbstract.toString(q));
    for (int i = 0; i < n; ++i) {
      assertEquals(n, q[i].length);
      Exam.assertDistribution(q[i]);
    }
    assertTrue(q[0][1] > q[0][2]); //""->A > ""->AA
    assertTrue(q[2][1] > q[2][0]); //AA->A > AA -> ""
  }

  public void testDiploid() throws InvalidParamsException, IOException {
    final DescriptionCommon desc = new DescriptionCommon("", "A", "AA");
    final Hypotheses<DescriptionCommon> hyp = new MockHypotheses<>(desc, SimplePossibility.SINGLETON, false, new double[] {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 0);
    //System.err.println(hyp.haploid());
    //System.err.println(hyp);
    final int n = hyp.size();
    final double[][] initialPriors = {
      {0.5, 0.3, 0.2},
      {0.25, 0.5, 0.25},
      {0.5 / 3.0, 1.0 / 3.0, 0.5},
    };
    final double[][] q = SomaticPriors.makeQ(0.3, 0.0, hyp, initialPriors);
    assertEquals(n, q.length);
    //System.err.println(IntegralAbstract.toString(q));
    for (int i = 0; i < n; ++i) {
      assertEquals(n, q[i].length);
      Exam.assertDistribution(q[i]);
    }
    //see spreadsheet
    final double[][] exp = {
      {0.723, 0.008,  0.004,  0.153,  0.011,  0.102},
      {0.006, 0.723,  0.006,  0.128,  0.128,  0.011},
      {0.003, 0.010,  0.723,  0.010,  0.170,  0.085},
      {0.064, 0.077,  0.005,  0.729,  0.058,  0.068},
      {0.004, 0.085,  0.064,  0.050,  0.730,  0.068},
      {0.043, 0.009,  0.051,  0.090,  0.083,  0.726},
    };
    checkEquals(q, exp);
  }

  private void checkEquals(final double[][] act, final double[][] exp) {
    final int len = exp.length;
    assertEquals(len, act.length);
    for (int i = 0; i < len; ++i) {
      final int leni = exp[i].length;
      assertEquals(leni, act[i].length);
      for (int j = 0; j < leni; ++j) {
        assertEquals(exp[i][j], act[i][j], 1e-3);
      }
    }
  }
}
