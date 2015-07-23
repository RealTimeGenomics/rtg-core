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

package com.rtg.variant.bayes;


import static com.rtg.util.StringUtils.LS;
import static com.rtg.util.StringUtils.TAB;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

import com.rtg.util.MathUtils;
import com.rtg.util.Utils;
import com.rtg.variant.StaticThreshold;
import com.rtg.variant.Variant;
import com.rtg.variant.VariantOutputLevel;
import com.rtg.variant.VariantParams;
import com.rtg.variant.VariantSample;
import com.rtg.variant.bayes.multisample.HaploidDiploidHypotheses;
import com.rtg.variant.bayes.multisample.singleton.SingletonCaller;
import com.rtg.variant.bayes.snp.DescriptionCommon;
import com.rtg.variant.bayes.snp.DescriptionSnp;
import com.rtg.variant.bayes.snp.EvidenceQ;
import com.rtg.variant.bayes.snp.HypothesesNone;
import com.rtg.variant.bayes.snp.HypothesesPrior;
import com.rtg.variant.bayes.snp.StatisticsSnp;
import com.rtg.variant.format.VariantOutputVcfFormatter;
import com.rtg.variant.util.VariantUtils;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;
import com.rtg.variant.util.arithmetic.SimplePossibility;

import junit.framework.TestCase;

/**
 */
public class ModelTest extends TestCase {

  @SuppressWarnings("unchecked")
  public static Variant makeCalls(final ModelInterface<?> mod, String refName, int start, int end, byte[] ref, final VariantParams outputParams) {
    final SingletonCaller sc = new SingletonCaller(outputParams);
    final Hypotheses<?> hyp = mod.hypotheses();
    final HaploidDiploidHypotheses<HypothesesPrior<Description>> hyp2 = new HaploidDiploidHypotheses<>(HypothesesNone.SINGLETON,
                                                                                                               (HypothesesPrior<Description>) (hyp.haploid() ? hyp : null),
                                                                                                               (HypothesesPrior<Description>) (hyp.haploid() ? null : hyp));
    final List<ModelInterface<?>> mods = new ArrayList<>();
    mods.add(mod);
    return sc.makeCall(refName, start, end, ref, mods, hyp2);
  }


  public void testAmbiguity() {
    for (int i = 0; i <= Model.AMBIGUITY_PHRED; i++) {
      assertTrue(VariantUtils.phredToProb(i) >= Model.AMBIGUITY_THRESHOLD);
    }
    assertEquals(VariantUtils.phredToProb(Model.AMBIGUITY_PHRED), Model.AMBIGUITY_THRESHOLD);
    assertTrue(VariantUtils.phredToProb(Model.AMBIGUITY_PHRED + 1) < Model.AMBIGUITY_THRESHOLD);
  }

  public void testPosterior0() {
    final double[] priors = {0.1, 0.4, 0.35, 0.15};
    final PossibilityArithmetic arith = SimplePossibility.SINGLETON;
    final MockHypotheses<DescriptionCommon> hypotheses = new MockHypotheses<DescriptionCommon>(DescriptionSnp.SINGLETON, arith, true, priors, 0);
    final Model<?> mo = new MockModel<>(hypotheses, new StatisticsSnp(hypotheses.description()), null);
    assertTrue(mo.hypotheses() == hypotheses);
    assertEquals(4, mo.size());

    for (int i = 0; i < 4; i++) {
      assertEquals(0.0, mo.posteriorLn0(i));
    }
    final VariantParams params0 = VariantParams.builder().create();
    final Variant call0 = makeCalls(mo, null, 42, 43, null, params0);
    assertNull(call0);
  }

  //no increment force the posteriors
  public void testPosterior1() {
    final double[] priors = {0.1, 0.4, 0.35, 0.15};
    final PossibilityArithmetic arith = SimplePossibility.SINGLETON;
    final MockHypotheses<DescriptionCommon> hypotheses = new MockHypotheses<DescriptionCommon>(DescriptionSnp.SINGLETON, arith, true, priors, 0);
    final double[] post = {0.1, 0.6, 0.25, 0.05};
    final Model<?> mo = new MockModel<>(hypotheses, new StatisticsSnp(hypotheses.description()), post);
    assertTrue(mo.hypotheses() == hypotheses);
    assertEquals(4, mo.size());

    for (int i = 0; i < 4; i++) {
      assertEquals(Math.log(post[i]), mo.posteriorLn0(i));
    }

    final VariantParams params = VariantParams.builder().callLevel(VariantOutputLevel.ALL).create();
    final Variant call = makeCalls(mo, "foo", 42, 43, new byte[] {}, params);
    assertNotNull(call);
    assertTrue(call.isInteresting());
    assertEquals(3.5, call.getNonIdentityPosterior(), 0.1);
    assertEquals("", call.getSample(0).getStatisticsString());
    assertFalse(call.isIndel());
    final String exp = ""
        + "Model" + LS
        + "A  -2.303" + LS
        + "C  -0.511" + LS
        + "G  -1.386" + LS
        + "T  -2.996" + LS
        ;
    assertEquals(exp, mo.toString());
  }

  //simple case - single increment
  public void testIncrement1a() {
    final double[] priors = {0.1, 0.4, 0.35, 0.15};
    final PossibilityArithmetic arith = SimplePossibility.SINGLETON;
    final MockHypotheses<DescriptionCommon> hypotheses = new MockHypotheses<DescriptionCommon>(DescriptionSnp.SINGLETON, arith, true, priors, 0);
    final ModelInterface<DescriptionCommon> mo = new Model<>(hypotheses, new StatisticsSnp(hypotheses.description()));

    final double[] prob = {0.175, 0.25, 0.5, 0.075};
    final EvidenceInterface di = new MockEvidence(DescriptionSnp.SINGLETON, 0.0, prob, 0);
    mo.increment(di);

    for (int i = 0; i < 4; i++) {
      assertEquals(Math.log(prob[i]), mo.posteriorLn0(i));
    }

    final VariantParams params = VariantParams.builder().callLevel(VariantOutputLevel.ALL).create();
    final Variant call = makeCalls(mo, "foo", 42, 43, new byte[] {}, params);
    assertNotNull(call);
    assertTrue(call.isInteresting());
    assertEquals(2.795, call.getNonIdentityPosterior(), 0.001);
    assertFalse(call.isIndel());
    //    assertTrue(call.isAllowOutput());

    //see numbers above
    final VariantSample vs = call.getSample(0);
    assertNotNull(vs.getName());
    assertEquals(0.307, vs.getPosterior(), 0.001);
    assertEquals("G", vs.getName());
    assertFalse(vs.isIdentity());
    assertEquals(VariantSample.DeNovoStatus.UNSPECIFIED, vs.isDeNovo());

    assertEquals(2.795, call.getNonIdentityPosterior(), 0.001);
  }

  //has to very carefully set the reference not to be 0
  public void testNonIdentityPosterior() {
    final double[] priors = {0.1, 0.4, 0.35, 0.15};
    final PossibilityArithmetic arith = SimplePossibility.SINGLETON;
    final MockHypotheses<DescriptionCommon> hypotheses = new MockHypotheses<DescriptionCommon>(DescriptionSnp.SINGLETON, arith, true, priors, 0) {
      @Override
      public int reference() {
        return 1;
      }
    };
    final ModelInterface<DescriptionCommon> mo = new Model<>(hypotheses, new StatisticsSnp(hypotheses.description()));

    final double[] prob = {0.175, 0.25, 0.5, 0.075};
    final EvidenceInterface di = new MockEvidence(DescriptionSnp.SINGLETON, 0.0, prob, 2);
    mo.increment(di);

    for (int i = 0; i < 4; i++) {
      assertEquals(Math.log(prob[i]), mo.posteriorLn0(i));
    }

    final VariantParams params = VariantParams.builder().callLevel(VariantOutputLevel.ALL).create();
    final Variant call = makeCalls(mo, "foo", 42, 43, new byte[] {}, params);
    assertNotNull(call);
    assertTrue(call.isInteresting());
    assertEquals(0.712, call.getNonIdentityPosterior(), 0.001);
    assertFalse(call.isIndel());
    //    assertTrue(call.isAllowOutput());

    //see numbers above
    final VariantSample vs = call.getSample(0);
    assertNotNull(vs.getName());
    assertEquals(0.307, vs.getPosterior(), 0.001);
    assertEquals("G", vs.getName());
    assertFalse(vs.isIdentity());
  }

  //simple single non-ref count
  public void testIncrement1b() {
    final double[] priors = {0.1, 0.4, 0.35, 0.15};
    final PossibilityArithmetic arith = SimplePossibility.SINGLETON;
    final MockHypotheses<DescriptionCommon> hypotheses = new MockHypotheses<DescriptionCommon>(DescriptionSnp.SINGLETON, arith, true, priors, 0);
    final ModelInterface<DescriptionCommon> mo = new Model<>(hypotheses, new StatisticsSnp(hypotheses.description()));

    final double[] prob = {0.175, 0.25, 0.5, 0.075};
    final EvidenceInterface di = new MockEvidence(DescriptionSnp.SINGLETON, 0.0, prob, 1);
    mo.increment(di);

    for (int i = 0; i < 4; i++) {
      assertEquals(Math.log(prob[i]), mo.posteriorLn0(i));
    }

    final VariantParams params = VariantParams.builder().callLevel(VariantOutputLevel.ALL).create();
    final Variant call = makeCalls(mo, "foo", 42, 43, new byte[] {}, params);
    assertNotNull(call);
    assertTrue(call.isInteresting());
    assertEquals(2.795, call.getNonIdentityPosterior(), 0.001);
    assertFalse(call.isIndel());
    //    assertTrue(call.isAllowOutput());

    //see numbers above
    final VariantSample vs = call.getSample(0);
    assertNotNull(vs.getName());
    assertEquals(0.307, vs.getPosterior(), 0.001);
    assertEquals("G", vs.getName());
    assertFalse(vs.isIdentity());
  }

  public void testIncrement2() {
    final double[] priors = {0.1, 0.4, 0.35, 0.15};
    final PossibilityArithmetic arith = SimplePossibility.SINGLETON;
    final MockHypotheses<DescriptionCommon> hypotheses = new MockHypotheses<DescriptionCommon>(DescriptionSnp.SINGLETON, arith, true, priors, 0);
    final ModelInterface<DescriptionCommon> mo = new Model<>(hypotheses, new StatisticsSnp(hypotheses.description()));

    final double[] prob = {0.175, 0.25, 0.5, 0.075};
    final EvidenceInterface di = new MockEvidence(DescriptionSnp.SINGLETON, 0.0, prob, 0);
    mo.increment(di);

    for (int i = 0; i < 4; i++) {
      assertEquals(Math.log(prob[i]), mo.posteriorLn0(i));
    }

    final VariantParams params = VariantParams.builder().callLevel(VariantOutputLevel.INTERESTING).create();
    final Variant call = makeCalls(mo, "foo", 42, 43, new byte[] {}, params);
    assertNull(call);
  }

  //increment twice set trigger - suppress complex calls
  public void testIncrementTriggerNoComplex() {
    final double[] priors = {0.1, 0.4, 0.35, 0.15};
    final PossibilityArithmetic arith = SimplePossibility.SINGLETON;
    final MockHypotheses<DescriptionCommon> hypotheses = new MockHypotheses<DescriptionCommon>(DescriptionSnp.SINGLETON, arith, true, priors, 0);
    final ModelInterface<DescriptionCommon> mo = new Model<>(hypotheses, new StatisticsSnp(hypotheses.description()));

    final double[] prob = {0.175, 0.25, 0.5, 0.075};
    final EvidenceInterface di = new MockEvidence(DescriptionSnp.SINGLETON, 0.0, prob, 1);
    mo.increment(di);
    mo.increment(di);

    for (int i = 0; i < 4; i++) {
      assertEquals(Math.log(prob[i] * prob[i]), mo.posteriorLn0(i));
    }

    final VariantParams params = VariantParams.builder().callLevel(VariantOutputLevel.ALL).noComplexCalls(true).create();
    final Variant call = makeCalls(mo, "foo", 42, 43, new byte[] {}, params);
    assertNotNull(call);
    assertTrue(call.isInteresting());
    assertEquals(3.612, call.getNonIdentityPosterior(), 0.001);
    assertFalse(call.isIndel());
    //    assertTrue(call.isAllowOutput());
    final VariantSample vs = call.getSample(0);
    assertNotNull(vs.getName());
    assertEquals(1.108, vs.getPosterior(), 0.001);
    assertEquals("G", vs.getName());
    assertFalse(vs.isIdentity());
  }

  public void testIncrementAmbiguous() {
    final double[] priors = {0.1, 0.4, 0.35, 0.15};
    final PossibilityArithmetic arith = SimplePossibility.SINGLETON;
    final MockHypotheses<DescriptionCommon> hypotheses = new MockHypotheses<DescriptionCommon>(DescriptionSnp.SINGLETON, arith, true, priors, 0);
    final ModelInterface<DescriptionCommon> mo = new Model<>(hypotheses, new StatisticsSnp(hypotheses.description()));

    final double[] prob = {0.175, 0.25, 0.5, 0.075};
    final double r = 0.49;
    mo.increment(new MockEvidence(DescriptionSnp.SINGLETON, r, prob, 1));
    mo.increment(new MockEvidence(DescriptionSnp.SINGLETON, 0.51, prob, 1));

    for (int i = 0; i < 4; i++) {
      final double p = prob[i] * (1.0 - r) + r * 0.25;
      assertEquals(Math.log(p), mo.posteriorLn0(i));
    }

    final VariantParams params = VariantParams.builder().callLevel(VariantOutputLevel.ALL).noComplexCalls(true).maxAmbiguity(0.5).create();
    final Variant call = makeCalls(mo, "foo", 42, 43, new byte[] {}, params);
    assertNotNull(call);
    //assertFalse(call.isInteresting());
    assertEquals(2.493, call.getNonIdentityPosterior(), 0.001);
    assertFalse(call.isIndel());
    //    assertTrue(call.isAllowOutput());
    final VariantSample vs = call.getSample(0);
    assertNotNull(vs.getName());
    assertEquals(-0.094, vs.getPosterior(), 0.001);
    assertEquals("G", vs.getName());
    assertFalse(vs.isIdentity());
  }

  public void testIncrementBoring() {
    final double[] priors = {0.4, 0.1, 0.35, 0.15};
    final PossibilityArithmetic arith = SimplePossibility.SINGLETON;
    final MockHypotheses<DescriptionCommon> hypotheses = new MockHypotheses<DescriptionCommon>(DescriptionSnp.SINGLETON, arith, true, priors, 0);
    final ModelInterface<DescriptionCommon> mo = new Model<>(hypotheses, new StatisticsSnp(hypotheses.description()));

    final double[] prob = {0.5, 0.175, 0.25, 0.075};
    final EvidenceInterface di = new MockEvidence(DescriptionSnp.SINGLETON, 0.0, prob, 0);
    mo.increment(di);
    mo.increment(di);
    final EvidenceInterface di1 = new MockEvidence(DescriptionSnp.SINGLETON, 0.0, prob, 1);
    mo.increment(di1);

    for (int i = 0; i < 4; i++) {
      final double cube = prob[i] * prob[i] * prob[i];
      assertEquals(Math.log(cube), mo.posteriorLn0(i));
    }

    final VariantParams params = VariantParams.builder().callLevel(VariantOutputLevel.INTERESTING).noComplexCalls(true).interestingThreshold(0.0).create();
    final Variant call = makeCalls(mo, "foo", 42, 43, new byte[] {}, params);
    assertNull(call);
  }

  //call equal to ref
  public void testIncrement3() {
    final double[] priors = {0.4, 0.1, 0.35, 0.15};
    final PossibilityArithmetic arith = SimplePossibility.SINGLETON;
    final MockHypotheses<DescriptionCommon> hypotheses = new MockHypotheses<DescriptionCommon>(DescriptionSnp.SINGLETON, arith, true, priors, 0);
    final ModelInterface<DescriptionCommon> mo = new Model<>(hypotheses, new StatisticsSnp(hypotheses.description()));

    final double[] prob = {0.5, 0.175, 0.25, 0.075};
    final EvidenceInterface di = new MockEvidence(DescriptionSnp.SINGLETON, 0.0, prob, 0);
    mo.increment(di);

    for (int i = 0; i < 4; i++) {
      assertEquals(Math.log(prob[i]), mo.posteriorLn0(i));
    }

    final VariantParams params = VariantParams.builder().callLevel(VariantOutputLevel.ALL).create();
    final Variant call = makeCalls(mo, "foo", 42, 43, new byte[] {}, params);
    assertNotNull(call);
    assertTrue(call.isInteresting());
    assertNotNull(call.getNonIdentityPosterior());
    assertFalse(call.isIndel());
    //    assertTrue(call.isAllowOutput());

    final VariantSample vs = call.getSample(0);
    assertNotNull(vs.getName());
    assertEquals("A", vs.getName());
    assertTrue(vs.isIdentity());
  }

  private static final DescriptionCommon MOCK_DESCRIPTION = new DescriptionCommon("X", "Y");


  //call not haploid
  public void testIncrement4() {
    final double[] priors = {0.1, 0.15, 0.75};
    final PossibilityArithmetic arith = SimplePossibility.SINGLETON;
    final MockHypotheses<DescriptionCommon> hypotheses = new MockHypotheses<>(MOCK_DESCRIPTION, arith, false, priors, 0);
    final ModelInterface<DescriptionCommon> mo = new Model<>(hypotheses, new StatisticsSnp(hypotheses.description()));

    final double[] prob = {0.8, 0.2};
    final EvidenceInterface di = new MockEvidence(DescriptionSnp.SINGLETON, 0.0, prob, 0);
    mo.increment(di);

    assertEquals(Math.log(prob[0]), mo.posteriorLn0(0));

    assertEquals(Math.log(prob[1]), mo.posteriorLn0(1));

    final double pr2 = (prob[0] + prob[1]) / 2.0;
    assertEquals(Math.log(pr2), mo.posteriorLn0(2));

    final VariantParams params = VariantParams.builder().callLevel(VariantOutputLevel.ALL).create();
    final Variant call = makeCalls(mo, "foo", 42, 43, new byte[] {}, params);
    assertNotNull(call);
    assertTrue(call.isInteresting());
    assertNotNull(call.getNonIdentityPosterior());
    assertFalse(call.isIndel());
    //    assertTrue(call.isAllowOutput());

    final VariantSample vs = call.getSample(0);
    assertNotNull(vs.getName());
    assertEquals("X:Y", vs.getName());
    assertFalse(vs.isIdentity());
  }

  private static final DescriptionCommon MOCK_PE_DESCRIPTION = new DescriptionCommon("X", "Y");

  //make r and pE not be 0
  public void testIncrement5() {
    final double[] priors = {0.6, 0.4};
    final PossibilityArithmetic arith = SimplePossibility.SINGLETON;
    final MockHypotheses<DescriptionCommon> hypotheses = new MockHypotheses<>(MOCK_PE_DESCRIPTION, arith, true, priors, 0);
    final ModelInterface<DescriptionCommon> mo = new Model<>(hypotheses, new StatisticsSnp(hypotheses.description()));

    final double[] prob = {0.75, 0.25};
    final EvidenceInterface di = new MockEvidence(DescriptionSnp.SINGLETON, 0.15, prob, 0)  {
      @Override
      public double pe() {
        return 0.1;
      }
    };
    mo.increment(di);

    assertEquals(-0.4269, mo.posteriorLn0(0), 0.0001);

    assertEquals(-1.4806, mo.posteriorLn0(1), 0.0001);
  }

  // COPIED FROM BayesianTest

  private static final Description OLD_DESCRIPTION = new DescriptionCommon("A", "BC", "BC", "BC");

  private static final class OldHypotheses extends HypothesesPrior<Description> {
    OldHypotheses(final double[] priors, int ref) {
      super(OLD_DESCRIPTION, SimplePossibility.SINGLETON, priors, true, ref);
    }
  }

  private static final VariantOutputVcfFormatter FORMATTER = new VariantOutputVcfFormatter();

  public void testNonRefCount0a() {
    final OldHypotheses hy = new OldHypotheses(new double[] {0.1, 0.2, 0.7}, 0);
    final MockModel<?> ba = new MockModel<>(hy, new StatisticsSnp(hy.description()), null);
    ba.integrity();

    assertEquals(0, ba.statistics().nonRefCount());
    assertNull(makeCalls(ba, null, 0, 1, new byte[] {}, VariantParams.builder().create()));
  }

  //check that nonref count is unchanged under increments by ref
  public void testNonRefCount1() {
    final OldHypotheses hy = new OldHypotheses(new double[] {0.1, 0.2, 0.7}, 0);
    final MockModel<?> ba = new MockModel<>(hy, new StatisticsSnp(hy.description()), null);
    ba.integrity();

    final EvidenceInterface di = new EvidenceQ(OLD_DESCRIPTION, 0, 0, 0, 0.1, 0.1, true, false, false, false);
    for (int i = 0; i < 4; i++) {
      ba.increment(di);
    }
    ba.integrity();

    assertEquals(0, ba.statistics().nonRefCount());
  }

  //add 4 non reference increments
  public void testNonRefCount2() {
    final OldHypotheses hy = new OldHypotheses(new double[] {0.1, 0.2, 0.7}, 0);
    final MockModel<?> ba = new MockModel<>(hy, new StatisticsSnp(hy.description()), null);
    ba.integrity();

    final EvidenceInterface di = new EvidenceQ(OLD_DESCRIPTION, 1, 0, 0, 0.1, 0.1, true, false, false, false);
    for (int i = 0; i < 4; i++) {
      ba.increment(di);
    }
    ba.integrity();

    assertEquals(4, ba.statistics().nonRefCount());
  }

  public void testCoverageThreshold() {
    final OldHypotheses hy = new OldHypotheses(new double[] {0.1, 0.1, 0.2, 0.7}, 0);
    final MockModel<?> ba = new MockModel<>(hy, new StatisticsSnp(hy.description()), null);
    ba.integrity();

    final EvidenceInterface di = new EvidenceQ(OLD_DESCRIPTION, 1, 0, 0, 0.1, 0.1, true, false, false, false);
    for (int i = 0; i < 4; i++) {
      ba.increment(di);
    }
    ba.integrity();

    final String prefix = "x";
    final VariantParams options = VariantParams.builder()
        .maxCoverageFilter(new StaticThreshold(4))
        .maxAmbiguity(0.1)
        .create();
    // within coverage threshold
    final Variant vo = makeCalls(ba, prefix, 0, 1, new byte[] {}, options);
    final String bostr = FORMATTER.formatCall(vo);
    //    assertTrue(vo.isAllowOutput());
    assertTrue(bostr, bostr.contains("PASS"));

    ba.increment(di);

    // outside coverage threshold
    final Variant vo2 = makeCalls(ba, prefix, 0, 1, new byte[] {}, options);
    assertTrue(vo2.isInteresting());
    final String bostr2 = FORMATTER.formatCall(vo2);
    assertTrue(bostr2, bostr2.contains("OC"));
    assertFalse(bostr2, bostr2.contains("PASS"));
  }

  public void testIHThreshold() {
    final OldHypotheses hy = new OldHypotheses(new double[] {0.1, 0.1, 0.2, 0.7}, 0);
    final MockModel<?> ba = new MockModel<>(hy, new StatisticsSnp(hy.description()), null);
    ba.integrity();

    final EvidenceInterface di = new EvidenceQ(OLD_DESCRIPTION, 1, 0, 0, 0.49, 0.1, true, false, false, false);
    for (int i = 0; i < 9; i++) {
      ba.increment(di);
    }
    final EvidenceInterface di1 = new EvidenceQ(OLD_DESCRIPTION, 1, 0, 0, 0.51, 0.1, true, false, false, false);
    ba.increment(di1);
    ba.integrity();

    final String prefix = "x";
    final VariantParams options = VariantParams.builder()
        .maxCoverageFilter(new StaticThreshold(100))
        .maxAmbiguity(0.1)
        .create();
    final VariantOutputVcfFormatter formatter = new VariantOutputVcfFormatter(options, "SAMPLE");

    // within IHRatio
    //System.err.println(IntegralAbstract.toString(ba));
    final Variant vo = makeCalls(ba, prefix, 0, 1, new byte[] {}, options);
    assertTrue(vo.isInteresting());
    final String bostr = formatter.formatCall(vo);
    assertTrue(bostr, bostr.contains("PASS"));
    ba.increment(di1);

    // outside IHRatio
    final Variant vo2 = makeCalls(ba, prefix, 0, 1, new byte[] {}, options);
    assertTrue(vo2.isInteresting());
    final String bostr2 = formatter.formatCall(vo2);
    assertTrue(bostr2, bostr2.contains("a10.0"));
    assertFalse(bostr2, bostr2.contains("PASS"));
  }


  public void testPosteriorThreshold2() {
    final OldHypotheses hy = new OldHypotheses(new double[] {0.01, 0.89, 0.05, 0.05}, 1);
    final MockModel<?> ba = new MockModel<>(hy, new StatisticsSnp(hy.description()), null);
    ba.integrity();

    final Variant vo = makeCalls(ba, "x", 0, 1, new byte[] {}, VariantParams.builder().create());
    assertNull(vo);
  }

  protected void checkQualityOutput(final double p, final String exp) {
    final ByteArrayOutputStream ba = new ByteArrayOutputStream();
    final PrintStream pr = new PrintStream(ba);
    final double ln = Math.log(p / (1.0 - p));
    final double ql = ln / MathUtils.LOG_10;
    final String fq = Utils.realFormat(ql, 1);
    pr.append("\t").append(fq);
    pr.close();
    assertEquals(exp, ba.toString());
  }

  public void testQualityOut() {
    checkQualityOutput(0.5, TAB + "0.0");
    checkQualityOutput(0.1, TAB + "-1.0");
    checkQualityOutput(0.2, TAB + "-0.6");
    checkQualityOutput(0.99, TAB + "2.0");
  }

  public void testAmbiguityShortcut() {
    assertFalse(new EvidenceQ(DescriptionSnp.SINGLETON, 1, 0, 0, 0.3, 0.7, true, false, false, false).mapError() >= Model.AMBIGUITY_THRESHOLD);
    assertTrue(new EvidenceQ(DescriptionSnp.SINGLETON, 1, 0, 0, 0.51, 0.8, true, false, false, false).mapError() >= Model.AMBIGUITY_THRESHOLD);
  }
}
