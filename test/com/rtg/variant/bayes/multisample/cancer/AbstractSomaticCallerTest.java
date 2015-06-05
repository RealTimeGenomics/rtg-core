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
import java.util.EnumSet;
import java.util.List;

import com.rtg.util.InvalidParamsException;
import com.rtg.variant.Variant;
import com.rtg.variant.VariantOutputLevel;
import com.rtg.variant.VariantParams;
import com.rtg.variant.VariantParamsBuilder;
import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.Evidence;
import com.rtg.variant.bayes.Hypotheses;
import com.rtg.variant.bayes.Model;
import com.rtg.variant.bayes.ModelInterface;
import com.rtg.variant.bayes.multisample.HaploidDiploidHypotheses;
import com.rtg.variant.bayes.snp.DescriptionCommon;
import com.rtg.variant.bayes.snp.DescriptionSnp;
import com.rtg.variant.bayes.snp.EvidenceQ;
import com.rtg.variant.bayes.snp.HypothesesNone;
import com.rtg.variant.bayes.snp.HypothesesPrior;
import com.rtg.variant.bayes.snp.StatisticsSnp;
import com.rtg.variant.dna.DNARange;
import com.rtg.variant.dna.DNARangeAT;
import com.rtg.variant.format.VariantOutputVcfFormatter;
import com.rtg.variant.format.VcfFormatField;
import com.rtg.variant.format.VcfInfoField;
import com.rtg.variant.util.arithmetic.SimplePossibility;

import junit.framework.TestCase;

/**
 * @param <D> description type
 */
public abstract class AbstractSomaticCallerTest<D extends Description> extends TestCase {

  protected static HypothesesPrior<Description> simpleHomoHyps(final double same, final int ref) {
    final DescriptionCommon desc = DescriptionSnp.SINGLETON;
    return new HypothesesPrior<Description>(desc, SimplePossibility.SINGLETON, true, ref) {

      @Override
      public double p(int hyp) {
        return hyp == ref ? same : (1.0 - same) / 3;
      }
      @Override
      public int reference() {
        return ref;
      }

    };
  }

  List<ModelInterface<Description>> getNormalModel() {
    final List<ModelInterface<Description>> models = new ArrayList<>();
    for (int ref = 0; ref < 4; ref++) {
      final Hypotheses<Description> hyps = simpleHomoHyps(0.99, ref);
      models.add(new Model<>(hyps, new StatisticsSnp(hyps.description())));
    }
    return models;
  }

  protected abstract List<ModelInterface<D>> getModel();

  protected abstract AbstractSomaticCaller getSomaticCaller(final double mutation, final Hypotheses<D> hypotheses, String normalName, String cancerName, VariantParams params);

  protected abstract Hypotheses<D> getCancerHypotheses(final double same, final int ref);

  protected VariantOutputVcfFormatter getFormatter() {
    final VariantOutputVcfFormatter formatter = new VariantOutputVcfFormatter("normal", "cancer");
    formatter.addExtraInfoFields(EnumSet.of(VcfInfoField.SOMATIC, VcfInfoField.LOH, VcfInfoField.NCS));
    formatter.addExtraFormatFields(EnumSet.of(VcfFormatField.SSC, VcfFormatField.SS));
    return formatter;
  }

  /**
   * Simulates a given number of reads, all containing the <code>readNt</code> nucleotide.
   * Uses and returns a Bayesian returned from <code>getBayes()</code>.
   *
   * @param numReads how many reads
   * @param readNt the nucleotide in the read at the current position.
   * @return the Bayesian after the reads.
   */
  protected final List<ModelInterface<D>> doReads(final int numReads, final int readNt) {
    final List<ModelInterface<D>> models = getModel();
    for (final ModelInterface<D> m : models) {
      for (int i = 0; i < numReads; i++) {
        final Evidence ev = new EvidenceQ(m.description(), readNt, 0, 0, 0.05, 0.05, true, false, false, false);
        m.increment(ev);
      }
    }
    return models;
  }

  protected List<ModelInterface<Description>> doNormalReads(final int numReads, final int readNt) {
    final List<ModelInterface<Description>> models = getNormalModel();
    for (final ModelInterface<Description> m : models) {
      for (int i = 0; i < numReads; i++) {
        final Evidence ev = new EvidenceQ(m.description(), readNt, 0, 0, 0.05, 0.05, true, false, false, false);
        m.increment(ev);
      }
    }
    return models;
  }

  // check that three A reads give us the expected input values.
  public void test3a() {
    assertEquals(""
        + "Hyp Post0 A C G T" + LS
        + "A 7.660e-01" + LS
        + "C 2.275e-05" + LS
        + "G 2.275e-05" + LS
        + "T 2.275e-05" + LS
        , dump(doNormalReads(3, DNARangeAT.A)));
  }

  // check that three C reads give us the expected input values.
  public void test3C() {
    assertEquals(""
        + "Hyp Post0 A C G T" + LS
        + "A 2.275e-05" + LS
        + "C 7.660e-01" + LS
        + "G 2.275e-05" + LS
        + "T 2.275e-05" + LS
        , dump(doNormalReads(3, DNARangeAT.C)));
  }

  // check that three G reads give us the expected input values.
  public void test3G() {
    assertEquals(""
        + "Hyp Post0 A C G T" + LS
        + "A 2.275e-05" + LS
        + "C 2.275e-05" + LS
        + "G 7.660e-01" + LS
        + "T 2.275e-05" + LS
        , dump(doNormalReads(3, DNARangeAT.G)));
  }

  protected static final String EXPECT_ALL_SAME = "chr1\t14\t.\tA\t.\t.\tPASS\tDP=6\tGT:DP:RE:AR:GQ:ABP:SBP:RPB:PUR:RS:AD:SS\t0:3:0.293:0.000:100:0.00:6.51:0.00:0.00:A,3,0.293:3\t0:3:0.293:0.000:75:0.00:6.51:0.00:0.00:A,3,0.293:3:0\n";

  public void testAllSame() throws InvalidParamsException, IOException {
    checkCancer(
      doNormalReads(3, DNARangeAT.A),
      doReads(3, DNARangeAT.A),
      EXPECT_ALL_SAME,
      "A",
      "C"
    );
  }

  protected static final String EXPECT_NORMAL_EQ_REF = "chr1\t14\t.\tA\tC\t.\tPASS\tSOMATIC=C;NCS=10.9;DP=6\tGT:DP:RE:AR:GQ:ABP:SBP:RPB:PUR:RS:AD:SSC:SS\t0:3:0.293:0.000:36:0.00:6.51:0.00:0.00:A,3,0.293:3,0\t1:3:0.293:0.000:11:0.00:6.51:0.00:0.00:C,3,0.293:0,3:1.0:2\n";

  public void testNormalEqualsRef() throws InvalidParamsException, IOException {
    checkCancer(
      doNormalReads(3, DNARangeAT.A),
      doReads(3, DNARangeAT.C),
      EXPECT_NORMAL_EQ_REF,
      "A",
      "C"
    );
  }

  protected static final String EXPECT_CANCER_EQ_REF = "chr1\t14\t.\tA\t.\t.\tPASS\tDP=6\tGT:DP:RE:AR:GQ:ABP:SBP:RPB:PUR:RS:AD:SS\t0:3:0.293:0.000:14:26.06:.:0.00:0.00:C,3,0.293:0\t0:3:0.293:0.000:25:0.00:6.51:0.00:0.00:A,3,0.293:3:0\n";

  public void testCancerEqualsRef() throws InvalidParamsException, IOException {
    checkCancer(
      doNormalReads(3, DNARangeAT.C),
      doReads(3, DNARangeAT.A),
      EXPECT_CANCER_EQ_REF,
      "A",
      "C"
    );
  }

  protected static final String EXPECT_CANCER_EQ_NORMAL = "chr1\t14\t.\tA\tC\t.\tPASS\tDP=6\tGT:DP:RE:AR:GQ:ABP:SBP:RPB:PUR:RS:AD:SS\t1:3:0.293:0.000:55:0.00:6.51:0.00:0.00:C,3,0.293:0,3\t1:3:0.293:0.000:65:0.00:6.51:0.00:0.00:C,3,0.293:0,3:1\n";

  public void testCancerEqualsNormal() throws InvalidParamsException, IOException {
    checkCancer(
      doNormalReads(3, DNARangeAT.C),
      doReads(3, DNARangeAT.C),
      EXPECT_CANCER_EQ_NORMAL,
      "A",
      "C"
    );
  }

  protected static final String EXPECT_ALL_DIFFERENT = "chr1\t14\t.\tA\tC,G\t.\tPASS\tSOMATIC=G;NCS=8.2;DP=6\tGT:DP:RE:AR:GQ:ABP:SBP:RPB:PUR:RS:AD:SSC:SS\t1:3:0.293:0.000:11:0.00:6.51:0.00:0.00:C,3,0.293:0,3,0\t2:3:0.293:0.000:11:0.00:6.51:0.00:0.00:G,3,0.293:0,0,3:0.7:2\n";

  public void testAllDifferent() throws InvalidParamsException, IOException {
    checkCancer(
      doNormalReads(3, DNARangeAT.C),
      doReads(3, DNARangeAT.G),
      EXPECT_ALL_DIFFERENT,
      "A",
      "C"
    );
  }

  private void checkCancer(List<ModelInterface<Description>> normal, List<ModelInterface<D>> cancer, String expect, String normalName, String cancerName) throws InvalidParamsException, IOException {
    final int refNt = DNARange.A;
    final int refCode = refNt - 1;
    final VariantParams params = new VariantParamsBuilder().callLevel(VariantOutputLevel.ALL).create();
    final Hypotheses<D> hypotheses = getCancerHypotheses(0.99, refCode);
    final AbstractSomaticCaller ccs = getSomaticCaller(0.001, hypotheses, normalName, cancerName, params);
    ccs.integrity();
    //System.out.println(new Posterior(ccs.mQ, A, normal.categories(), cancer.categories()));
    //    final MultivarianceCall out = new MultivarianceCall("chr1", 13, 14, VarianceCallType.SNP);
    //    out.setName("foo");
    final byte[] ref = new byte[14];
    ref[12] = DNARange.G;
    ref[13] = refNt;
    //debugCalls(outParams, 13, normal[refCode], cancer[refCode]);

    final List<ModelInterface<?>> models = new ArrayList<>();
    models.add(normal.get(refCode));
    models.add(cancer.get(refCode));
    final HaploidDiploidHypotheses<HypothesesPrior<Description>> hyp = new HaploidDiploidHypotheses<>(HypothesesNone.SINGLETON,
        (HypothesesPrior<Description>) normal.get(refCode).hypotheses(),
        null);
    final Variant v = ccs.makeCall("chr1", 13, 14, ref, models, hyp);
    final VariantOutputVcfFormatter formatter = getFormatter();
    assertEquals(expect, formatter.formatCall(v));
  }

  protected static <T extends Description> String dump(List<ModelInterface<T>> bayes) {
    final StringBuilder sb = new StringBuilder();
    sb.append("Hyp Post0 A C G T").append(LS);
    for (int i = 0; i < 4; i++) {
      sb.append(bayes.get(0).name(i)).append(" ").append(formatExp(Math.exp(bayes.get(0).posteriorLn0(i))));
      sb.append(LS);
    }
    return sb.toString();
  }

  protected static String formatExp(double value) {
    return String.format("%.3e", value);
  }
}
