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
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import com.rtg.reference.Ploidy;
import com.rtg.util.InvalidParamsException;
import com.rtg.variant.GenomePriorParams;
import com.rtg.variant.Variant;
import com.rtg.variant.VariantOutputLevel;
import com.rtg.variant.VariantParams;
import com.rtg.variant.VariantParamsBuilder;
import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.Evidence;
import com.rtg.variant.bayes.Hypotheses;
import com.rtg.variant.bayes.Model;
import com.rtg.variant.bayes.ModelInterface;
import com.rtg.variant.bayes.NoAlleleBalance;
import com.rtg.variant.bayes.complex.StatisticsComplex;
import com.rtg.variant.bayes.multisample.HaploidDiploidHypotheses;
import com.rtg.variant.bayes.snp.DescriptionCommon;
import com.rtg.variant.bayes.snp.EvidenceQ;
import com.rtg.variant.bayes.snp.HypothesesNone;
import com.rtg.variant.bayes.snp.HypothesesPrior;
import com.rtg.variant.bayes.snp.HypothesesSnp;
import com.rtg.variant.bayes.snp.StatisticsSnp;
import com.rtg.variant.dna.DNARange;
import com.rtg.variant.dna.DNARangeAT;
import com.rtg.variant.format.VariantOutputVcfFormatter;
import com.rtg.variant.util.arithmetic.LogApproximatePossibility;
import com.rtg.variant.util.arithmetic.SimplePossibility;

/**
 * This just contains a few extra tests that are specific to <code>PureSomaticCaller</code>.
 * They may be subsumed by the superclass tests.
 */
public class PureSomaticCallerTest extends AbstractSomaticCallerTest<Description> {

  @Override
  protected List<ModelInterface<Description>> getModel(Ploidy ploidy, double contamination, double same) {
    if (contamination > 0.0) {
      throw new IllegalArgumentException("Contamination should be 0 for pure somatic tests");
    }
    return getNormalModel(ploidy, 0.99);
  }

  @Override
  protected AbstractSomaticCaller getSomaticCaller(final Hypotheses<Description> hypotheses, VariantParams params, double phi, double psi) {
    return new PureSomaticCaller(new DefaultSomaticPriorsFactory<>(hypotheses, 0), new DefaultSomaticPriorsFactory<>(hypotheses, 0), params, phi, psi);
  }

  /** The result of 3 A reads, when reference is also A. */
  static final List<ModelInterface<Description>> EQUALS_REF_A =
      new PureSomaticCallerTest().doCancerReads(3, DNARangeAT.A);

  /** The result of 3 C reads, when reference is A. */
  static final List<ModelInterface<Description>> SEEN_3_C =
      new PureSomaticCallerTest().doCancerReads(3, DNARangeAT.C);

  /** The result of 3 G reads, when reference is A. */
  static final List<ModelInterface<Description>> SEEN_3_G =
      new PureSomaticCallerTest().doCancerReads(3, DNARangeAT.G);

  protected static final String EXPECT_IDENTICAL = "chr1\t14\t.\tG\tA\t.\tPASS\tDP=2\tGT:DP:RE:AR:GQ:ABP:SBP:RPB:AQ:PUR:RS:AD:SS\t1:1:0.020:0.000:21:0.00:2.17:0.00:0.000,17.011:0.00:A,1,0.020:0,1\t1:1:0.020:0.000:25:0.00:2.17:0.00:0.000,17.011:0.00:A,1,0.020:0,1:1\n";

  /**
   * Test that two identical SNP calls are not viewed as cancer.
   * @throws InvalidParamsException
   * @throws IOException
   */
  public void testIdentical() throws InvalidParamsException, IOException {
    final byte refNt = DNARange.G;
    final int refCode = refNt - 1;
    final HypothesesPrior<Description> hypotheses = simpleHyps(0.7, refCode, Ploidy.HAPLOID);
    final ModelInterface<?> model0 = new Model<>(hypotheses, new StatisticsSnp(hypotheses.description()), new NoAlleleBalance());
    final ModelInterface<?> model1 = new Model<>(hypotheses, new StatisticsSnp(hypotheses.description()), new NoAlleleBalance());
    final byte[] ref = new byte[14];
    ref[12] = DNARange.T;
    ref[13] = refNt;
    final Evidence ev = new EvidenceQ(hypotheses.description(), 0, 0, 0, 0.01, 0.01, true, false, false, false);
    model0.increment(ev);
    model1.increment(ev);
    final VariantParams params = new VariantParamsBuilder().callLevel(VariantOutputLevel.ALL).create();
    final AbstractSomaticCaller ccs = getSomaticCaller(simpleHyps(0.7, refCode, Ploidy.HAPLOID), params);
    ccs.integrity();
    final List<ModelInterface<?>> models = new ArrayList<>();
    models.add(model0);
    models.add(model1);
    for (ModelInterface<?> model : models) {
      model.freeze();
    }
    final Variant v = ccs.makeCall("chr1", 13, 14, ref, models, new HaploidDiploidHypotheses<>(HypothesesNone.SINGLETON, hypotheses, null));
    final VariantOutputVcfFormatter formatter = getFormatter();
    assertEquals(EXPECT_IDENTICAL, formatter.formatCall(v));
  }

  protected static final String EXPECT_CANCER1 = "chr1\t14\t.\tG\tA,C\t.\tPASS\tNCS=314.906;DP=30\tGT:DP:RE:AR:GQ:ABP:SBP:RPB:AQ:PUR:RS:AD:SSC:SS\t1:20:0.040:0.000:635:0.00:43.43:0.00:0.000,539.837,0.000:0.00:A,20,0.040:0,20,0\t2:10:0.020:0.000:314:0.00:21.71:0.00:0.000,0.000,269.919:0.00:C,10,0.020:0,0,10:31.4:2\n";

  /**
   * Test that a successful cancer call is made.
   * @throws InvalidParamsException
   * @throws IOException
   */
  public void testCancer1() throws InvalidParamsException, IOException {
    final byte refNt = DNARange.G;
    final int refCode = refNt - 1;
    final HypothesesPrior<Description> hypotheses = simpleHyps(0.7, refCode, Ploidy.HAPLOID);
    final ModelInterface<?> model0 = new Model<>(hypotheses, new StatisticsSnp(hypotheses.description()), new NoAlleleBalance());
    final ModelInterface<?> model1 = new Model<>(hypotheses, new StatisticsSnp(hypotheses.description()), new NoAlleleBalance());
    final byte[] ref = new byte[14];
    ref[12] = DNARange.T;
    ref[13] = refNt;

    final Evidence eva = new EvidenceQ(hypotheses.description(), 0, 0, 0, 0.001, 0.001, true, false, false, false);
    final Evidence evc = new EvidenceQ(hypotheses.description(), 1, 0, 0, 0.001, 0.001, true, false, false, false);

    for (int coverage = 0; coverage < 10; ++coverage) {
      model0.increment(eva);
      model0.increment(eva);
      model1.increment(evc);
    }
    final VariantParams params = new VariantParamsBuilder().callLevel(VariantOutputLevel.ALL).create();
    final AbstractSomaticCaller ccs = getSomaticCaller(hypotheses, params);
    ccs.integrity();
    final List<ModelInterface<?>> models = new ArrayList<>();
    models.add(model0);
    models.add(model1);
    for (ModelInterface<?> model : models) {
      model.freeze();
    }
    final Variant v = ccs.makeCall("chr1", 13, 14, ref, models, new HaploidDiploidHypotheses<>(HypothesesNone.SINGLETON, hypotheses, null));
    final VariantOutputVcfFormatter formatter = getFormatter();
    assertEquals(EXPECT_CANCER1, formatter.formatCall(v));
  }

  public void testQ() {
    final GenomePriorParams params = GenomePriorParams.builder().create();
    final VariantParams vParams = VariantParams.builder().genomePriors(params).create();
    final AbstractSomaticCaller ccs = getSomaticCaller(new HypothesesSnp(SimplePossibility.SINGLETON, params, false, 0), vParams);
    ccs.integrity();
    final String exp = ""
        + "length=10" + LS
        //   A     C     G     T   A:C   A:G   A:T   C:G   C:T   G:T
        + "0.490 0.010 0.010 0.010 0.140 0.020 0.020 0.140 0.020 0.140 " + LS
        + "0.010 0.490 0.010 0.010 0.140 0.140 0.020 0.020 0.140 0.020 " + LS
        + "0.010 0.010 0.490 0.010 0.020 0.140 0.140 0.140 0.020 0.020 " + LS
        + "0.010 0.010 0.010 0.490 0.020 0.020 0.140 0.020 0.140 0.140 " + LS
        + "0.070 0.070 0.010 0.010 0.500 0.080 0.020 0.080 0.080 0.080 " + LS
        + "0.010 0.070 0.070 0.010 0.080 0.500 0.080 0.080 0.080 0.020 " + LS
        + "0.010 0.010 0.070 0.070 0.020 0.080 0.500 0.080 0.080 0.080 " + LS
        + "0.070 0.010 0.070 0.010 0.080 0.080 0.080 0.500 0.020 0.080 " + LS
        + "0.010 0.070 0.010 0.070 0.080 0.080 0.080 0.020 0.500 0.080 " + LS
        + "0.070 0.010 0.010 0.070 0.080 0.020 0.080 0.080 0.080 0.500 " + LS
        ;
    assertEquals(exp, ccs.toString());
  }

  //this demonstrates the situation where we call 1/1 -> 1/2 (normal -> cancer) even when no evidence is provided
  //for the normal and contamination is zero
  public void testXRXPureCall() {
    final DescriptionCommon desc = new DescriptionCommon("ref", "A", "B"); //to resemble what evidence complex does
    final double[] priors = {0.0, -40.48, -33.38, -40.358, -69.708, -33.262};
    final Hypotheses<Description> hyp = new HypothesesPrior<>(desc, LogApproximatePossibility.SINGLETON, priors, false, 0);
    final Model<Description> modelNormal = new Model<>(hyp, new StatisticsComplex(desc, 1), new NoAlleleBalance());
    final Model<Description> modelCancer = new Model<>(hyp, new StatisticsComplex(desc, 1), new NoAlleleBalance());

    final ModelIncrementer<Description> cancerIncrementer = new ModelIncrementer<>(Collections.singletonList(modelCancer));

    final Variant var = getVariant(new ModelIncrementer<>(Collections.singletonList(modelNormal)).freeze(),
      cancerIncrementer
        .doReadsCustom(11, desc, 0, new double[] {0.333, 0.333, 0.333}, 0.3333, 0.6666, 1e-7)
        .doReadsCustom(2, desc, 2, new double[] {0.499, 1e-10, 0.499}, 0.4986, 0.51, 1e-7)
        .doReadsCustom(14, desc, 1, new double[] {1e-10, 0.995, 1e-10}, 0.002216, 0.998, 1e-7)
        .doReadsCustom(9, desc, 2, new double[] {1e-10, 1e-10, 0.995}, 0.00149, 0.998, 1e-7)
        .freeze(), getDefaultParams(), 0.01, 0.01);

    assertNotNull(var.getNormalCancerScore());
    assertEquals("B:B", var.getSample(0).getName());
    assertEquals("A:B", var.getSample(1).getName());
  }
}
