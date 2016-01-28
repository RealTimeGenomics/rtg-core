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

package com.rtg.variant.bayes.multisample.forwardbackward;


import com.rtg.variant.GenomePriorParams;
import com.rtg.variant.bayes.Code;
import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.Factor;
import com.rtg.variant.bayes.Hypotheses;
import com.rtg.variant.bayes.MockHypotheses;
import com.rtg.variant.bayes.Model;
import com.rtg.variant.bayes.ModelInterface;
import com.rtg.variant.bayes.NoAlleleBalance;
import com.rtg.variant.bayes.multisample.family.MendelianAlleleProbability;
import com.rtg.variant.bayes.multisample.family.MendelianAlleleProbabilityDiploid;
import com.rtg.variant.bayes.snp.DescriptionCommon;
import com.rtg.variant.bayes.snp.DescriptionNone;
import com.rtg.variant.bayes.snp.HypothesesNone;
import com.rtg.variant.bayes.snp.HypothesesPrior;
import com.rtg.variant.bayes.snp.HypothesesSnp;
import com.rtg.variant.bayes.snp.StatisticsSnp;
import com.rtg.variant.util.arithmetic.LogPossibility;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;
import com.rtg.variant.util.arithmetic.SimplePossibility;

import junit.framework.TestCase;

/**
 */
public class CommonFormulasTest extends TestCase {

  //  private void dumpM(final Code code, final MendelianAlleleProbability m) {
  //    final int range = 3;
  //    for (int k = 0; k < range; k++) {
  //      System.out.print(k);
  //      for (int h = 0; h < range; h++) {
  //        if (h != 0) {
  //          System.out.print(" ");
  //        }
  //        for (int j = 0; j < range; j++) {
  //          final double t = Math.exp(m.probabilityLn(code, j, k, h));
  //          System.out.format(" %1.2f", t);
  //        }
  //        System.out.println();
  //      }
  //      System.out.println();
  //    }
  //  }

  private static double[][] dm(final Factor<?> d, final MendelianAlleleProbability m) {
    final Hypotheses<?> hyp = d.hypotheses();
    final PossibilityArithmetic arith = d.arithmetic();
    final Code code = hyp.code();
    final int range = code.size();
    final double[][] dm = new double[range][range];
    for (int j = 0; j < range; j++) {
      for (int k = 0; k < range; k++) {
        double v = arith.zero();
        for (int h = 0; h < range; h++) {
          final double t = arith.multiply(d.p(h), arith.ln2Poss(m.probabilityLn(code, j, k, h)));
          v = arith.add(v, t);
        }
        dm[j][k] = v;
      }
    }
    return dm;
  }

  public void testForwardA() {
    final MendelianAlleleProbability m = MendelianAlleleProbabilityDiploid.SINGLETON;
    final Description desc = new DescriptionCommon("A", "B");
    final PossibilityArithmetic arith = SimplePossibility.SINGLETON;
    final Hypotheses<?> hyp = new MockHypotheses<>(desc, arith, false, new double[] {0.2, 0.3, 0.5}, 0);
    //dumpM(hyp.code(), m);
    final Factor<?> eu = new MutableFactor<>(hyp, arith, new double[] {0.1, 0.2, 0.3});
    final Factor<?> ev = new MutableFactor<>(hyp, arith, new double[] {0.15, 0.25, 0.35});

    final Factor<?> d0 = new MutableFactor<>(hyp, arith, new double[] {0.22, 0.33, 0.44});
    final Factor<?> d1 = new MutableFactor<>(hyp, arith, new double[] {0.11, 0.21, 0.31});
    final double[][][] c = {dm(d0, m), dm(d1, m)};
    final Factor<?> a = CommonFormulas.forwardA(eu, ev, hyp, 1, c, m);
    assertEquals(0.0259, a.p(0), 0.0001);
    assertEquals(0.0538, a.p(1), 0.0001);
    assertEquals(0.0841, a.p(2), 0.0001);
  }

  public void testDot() {
    final Description desc = new DescriptionCommon("A", "B");
    final PossibilityArithmetic arith = SimplePossibility.SINGLETON;
    final Hypotheses<?> hyp = new MockHypotheses<>(desc, arith, false, new double[] {0.2, 0.3, 0.5}, 0);
    //dumpM(hyp.code(), m);
    final Factor<?> a = new MutableFactor<>(hyp, arith, new double[] {0.1, 0.2, 0.3});
    final Factor<?> b = new MutableFactor<>(hyp, arith, new double[] {0.15, 0.25, 0.35});
    final Factor<?> c = CommonFormulas.dot(a, b);
    assertEquals(0.015, c.p(0));
    assertEquals(0.05, c.p(1));
    assertEquals(0.105, c.p(2));
  }

  public void testForwardE() {
    final Description desc = new DescriptionCommon("A", "B");
    final PossibilityArithmetic arith = SimplePossibility.SINGLETON;
    final Hypotheses<?> hyp = new MockHypotheses<>(desc, arith, false, new double[] {0.2, 0.3, 0.5}, 0);
    //dumpM(hyp.code(), m);
    final Factor<?> a = new MutableFactor<>(hyp, arith, new double[] {0.1, 0.2, 0.3});
    final Factor<?> b = new MutableFactor<>(hyp, arith, new double[] {0.15, 0.25, 0.35});
    //    final Factor<?> c = CommonFormulas.forwardE(a, b);
    //    final Factor<?> a = new MutableFactor<?>(hyp, new double[] {0.1, 0.2, 0.3});
    //    final Factor<?> b = new MutableFactor<?>(hyp, new double[] {0.15, 0.25, 0.35});
    final Factor<?> c = CommonFormulas.forwardEMonogamous(a, b);
    assertEquals(0.015, c.p(0));
    assertEquals(0.05, c.p(1));
    assertEquals(0.105, c.p(2));
  }

  public void testBackwardC() {
    final MendelianAlleleProbability m = MendelianAlleleProbabilityDiploid.SINGLETON;
    final Description desc = new DescriptionCommon("A", "B");
    final PossibilityArithmetic arith = SimplePossibility.SINGLETON;
    final Hypotheses<?> hyp = new MockHypotheses<>(desc, arith, false, new double[] {0.2, 0.3, 0.5}, 0);
    final Factor<?> a = new MutableFactor<>(hyp, arith, new double[] {0.1, 0.2, 0.3});
    final Factor<?> s = new MutableFactor<>(hyp, arith, new double[] {0.15, 0.25, 0.35});
    final Factor<?> b = new MutableFactor<>(hyp, arith, new double[] {0.22, 0.33, 0.44});

    final Factor<?> d = CommonFormulas.dot(s, b);
    assertEquals(0.15 * 0.22, d.p(0));
    assertEquals(0.25 * 0.33, d.p(1));
    assertEquals(0.35 * 0.44, d.p(2));
    final Factor<?> q = CommonFormulas.dot(a, d);
    assertEquals(0.1 * 0.15 * 0.22, q.p(0), 1e-8);
    assertEquals(0.2 * 0.25 * 0.33, q.p(1), 1e-8);
    assertEquals(0.3 * 0.35 * 0.44, q.p(2), 1e-8);

    final double[][] c = CommonFormulas.backwardC(hyp.size(), hyp.size(), CommonFormulas.maxCode(hyp, hyp),  d, m);
    //System.err.println(IntegralAbstract.toString(c));
    assertEquals(0.0330, c[0][0], 1e-4);
    assertEquals(0.1540, c[1][0], 1e-4);
    assertEquals(0.0935, c[2][0], 1e-4);
    assertEquals(0.1540, c[0][1], 1e-4);
    assertEquals(0.0825, c[1][1], 1e-4);
    assertEquals(0.1183, c[2][1], 1e-4);
    assertEquals(0.0935, c[0][2], 1e-4);
    assertEquals(0.1183, c[1][2], 1e-4);
    assertEquals(0.1059, c[2][2], 1e-4);
  }

  public void testBackwardBOneChild() {
    final MendelianAlleleProbability m = MendelianAlleleProbabilityDiploid.SINGLETON;
    final Description desc = new DescriptionCommon("A", "B");
    final PossibilityArithmetic arith = SimplePossibility.SINGLETON;
    final Hypotheses<?> hyp = new MockHypotheses<>(desc, arith, false, new double[] {0.2, 0.3, 0.5}, 0);
    final Factor<?> s = new MutableFactor<>(hyp, arith, new double[] {0.15, 0.25, 0.35});
    final Factor<?> b = new MutableFactor<>(hyp, arith, new double[] {0.22, 0.33, 0.44});
    final Factor<?> d = CommonFormulas.dot(s, b);

    final double[][] c0 = CommonFormulas.backwardC(hyp.size(), hyp.size(), CommonFormulas.maxCode(hyp, hyp), d, m);
    final double[][][] c = {c0};
    final Factor<?> ev = new MutableFactor<>(hyp, arith, new double[] {0.12, 0.13, 0.14});
    final Factor<?> backwardB = CommonFormulas.backwardB(s, ev, c, true);
    assertEquals(0.0371, backwardB.p(0), 1e-4);
    assertEquals(0.0458, backwardB.p(1), 1e-4);
    assertEquals(0.0414, backwardB.p(2), 1e-4);
  }

  public void testBackwardBTwoChildren() {
    final MendelianAlleleProbability m = MendelianAlleleProbabilityDiploid.SINGLETON;
    final Description desc = new DescriptionCommon("A", "B");
    final PossibilityArithmetic arith = SimplePossibility.SINGLETON;
    final Hypotheses<?> hyp = new MockHypotheses<>(desc, arith, false, new double[] {0.2, 0.3, 0.5}, 0);
    final Factor<?> s = new MutableFactor<>(hyp, arith, new double[] {0.15, 0.25, 0.35});
    final Factor<?> b = new MutableFactor<>(hyp, arith, new double[] {0.22, 0.33, 0.44});
    final Factor<?> d = CommonFormulas.dot(s, b);

    final double[][] c0 = CommonFormulas.backwardC(hyp.size(), hyp.size(), CommonFormulas.maxCode(hyp, hyp), d, m);
    final double[][][] c = {c0, c0};
    final Factor<?> ev = new MutableFactor<>(hyp, arith, new double[] {0.12, 0.13, 0.14});
    final Factor<?> backwardB = CommonFormulas.backwardB(s, ev, c, true);
    assertEquals(0.004438, backwardB.p(0), 1e-6);
    assertEquals(0.005688, backwardB.p(1), 1e-6);
    assertEquals(0.004436, backwardB.p(2), 1e-6);
  }

  public void testMaxCode() {
    final Description desc = new DescriptionCommon("A", "B");
    final double[] priors = {0.2, 0.3, 0.5};
    final Hypotheses<?> none = new HypothesesNone<>(DescriptionNone.SINGLETON, SimplePossibility.SINGLETON, 0);
    final Hypotheses<?> haploid = new MockHypotheses<>(desc, SimplePossibility.SINGLETON, true, priors, 0);
    final Hypotheses<?> diploid = new MockHypotheses<>(desc, SimplePossibility.SINGLETON, false, priors, 0);
    checkMaxCode(none, none, 0); //strictly this isn't needed as both parents cant be empty
    checkMaxCode(none, haploid, 2);
    checkMaxCode(haploid, haploid, 2);
    checkMaxCode(none, diploid, 3);
    checkMaxCode(haploid, diploid, 3);
    checkMaxCode(diploid, diploid, 3);
  }

  private void checkMaxCode(final Hypotheses<?> a, final Hypotheses<?> b, final int cnt) {
    assertEquals(cnt, CommonFormulas.maxCode(a, b).size());
    assertEquals(cnt, CommonFormulas.maxCode(b, a).size());
  }

  private static class MockModel extends Model<Description> {
    MockModel(final Hypotheses<Description> hyp, final double[] posteriors) {
      super(hyp, new StatisticsSnp(hyp.description()), new NoAlleleBalance());
      assert posteriors.length == mPosteriors.length;
      for (int i = 0; i < posteriors.length; i++) {
        mPosteriors[i] = arithmetic().prob2Poss(posteriors[i]);
      }
    }
  }

  public void testModelToVector() {
    final Description desc = new DescriptionCommon("A", "B");
    final Hypotheses<Description> hyp = new MockHypotheses<>(desc, SimplePossibility.SINGLETON, false, new double[] {0.0, 0.0, 0.0}/*not used*/, 0);
    final PossibilityArithmetic arith = hyp.arithmetic();
    final ModelInterface<?> model = new MockModel(hyp, new double[] {0.1, 0.2, 0.3});
    final Factor<?> ve = CommonFormulas.createMutableFactor(model);
    assertEquals(0.1, arith.poss2Prob(ve.p(0)), 1e-8);
    assertEquals(0.2, arith.poss2Prob(ve.p(1)), 1e-8);
    assertEquals(0.3, arith.poss2Prob(ve.p(2)), 1e-8);
  }

  public void testHypothesesToVector() {
    final Description desc = new DescriptionCommon("A", "B");
    final HypothesesPrior<?> hyp = new MockHypotheses<>(desc, SimplePossibility.SINGLETON, false, new double[] {0.1, 0.2, 0.3}, 0);
    final PossibilityArithmetic arith = hyp.arithmetic();
    final Factor<?> ve = CommonFormulas.createMutableFactor(hyp);
    assertEquals(0.1, arith.poss2Prob(ve.p(0)), 1e-8);
    assertEquals(0.2, arith.poss2Prob(ve.p(1)), 1e-8);
    assertEquals(0.3, arith.poss2Prob(ve.p(2)), 1e-8);
  }

  //Taken from HaploidSnpBayesianTest.testNewPriors()
  public void testPrior2VectorHaploid() {
    final GenomePriorParams params = GenomePriorParams.builder().create();
    final PossibilityArithmetic arith = LogPossibility.SINGLETON;
    final HypothesesSnp hyp = new HypothesesSnp(arith, params, true, 1);
    final Factor<?> hv = CommonFormulas.prior2Vector(params, hyp, arith);

    assertEquals(-9.257, hv.p(0), 0.001);
    assertEquals(-0.000, hv.p(1), 0.001);
    assertEquals(-9.172, hv.p(2), 0.001);
    assertEquals(-7.884, hv.p(3), 0.001);
  }

  //Taken from HaploidSnpBayesianTest.testNewPriors()
  public void testPrior2VectorDiploid() {
    final GenomePriorParams params = GenomePriorParams.builder().create();
    final PossibilityArithmetic arith = LogPossibility.SINGLETON;
    final HypothesesSnp hyp = new HypothesesSnp(arith, params, false, 1);
    final Factor<?> hv = CommonFormulas.prior2Vector(params, hyp, arith);

    assertEquals(-9.257, hv.p(0), 0.001);
    assertEquals(-0.002, hv.p(1), 0.001);
    assertEquals(-9.172, hv.p(2), 0.001);
    assertEquals(-7.885, hv.p(3), 0.001);

    assertEquals(-8.724,  hv.p(4), 0.001); //A:C
    assertEquals(-8.724,  hv.p(5), 0.001); //C:G
    assertEquals(-14.894, hv.p(6), 0.001); //G:T
    assertEquals(-14.894, hv.p(7), 0.001); //A:G
    assertEquals(-7.288,  hv.p(8), 0.001); //C:T
    assertEquals(-14.894, hv.p(9), 0.001); //A:T
  }

}
