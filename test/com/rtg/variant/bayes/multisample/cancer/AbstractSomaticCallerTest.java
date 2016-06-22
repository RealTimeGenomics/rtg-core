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

import com.rtg.reference.Ploidy;
import com.rtg.relation.GenomeRelationships;
import com.rtg.util.InvalidParamsException;
import com.rtg.variant.SomaticParamsBuilder;
import com.rtg.variant.Variant;
import com.rtg.variant.VariantOutputLevel;
import com.rtg.variant.VariantParams;
import com.rtg.variant.VariantParamsBuilder;
import com.rtg.variant.VariantSample;
import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.Hypotheses;
import com.rtg.variant.bayes.Model;
import com.rtg.variant.bayes.ModelInterface;
import com.rtg.variant.bayes.NoAlleleBalance;
import com.rtg.variant.bayes.multisample.HaploidDiploidHypotheses;
import com.rtg.variant.bayes.snp.DescriptionCommon;
import com.rtg.variant.bayes.snp.DescriptionSnp;
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

  protected static HypothesesPrior<Description> simpleHyps(final double same, final int ref, Ploidy ploidy) {
    final DescriptionCommon desc = DescriptionSnp.SINGLETON;
    return getDescriptionHypothesesPrior(same, ref, ploidy, desc);
  }

  protected static HypothesesPrior<Description> getDescriptionHypothesesPrior(final double same, final int ref, final Ploidy ploidy, final DescriptionCommon desc) {
    return new HypothesesPrior<Description>(desc, SimplePossibility.SINGLETON, ploidy == Ploidy.HAPLOID, ref) {

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

  ModelIncrementer<Description> getNormalIncremeter(Ploidy ploidy) {
    return getNormalIncremeter(ploidy, 0.99);
  }
  ModelIncrementer<Description> getNormalIncremeter(Ploidy ploidy, double same) {
    return new ModelIncrementer<>(getNormalModel(ploidy, same));
  }
  ModelIncrementer<D> getIncrementer(Ploidy ploidy, double contamination, double same) {
    return new ModelIncrementer<>(getModel(ploidy, contamination, same));
  }

  List<ModelInterface<Description>> getNormalModel(Ploidy ploidy, double same) {
    final List<ModelInterface<Description>> models = new ArrayList<>();
    for (int ref = 0; ref < 4; ref++) {
      final Hypotheses<Description> hyps = simpleHyps(same, ref, ploidy);
      models.add(new Model<>(hyps, new StatisticsSnp(hyps.description()), new NoAlleleBalance()));
    }
    return models;
  }

  protected abstract List<ModelInterface<D>> getModel(Ploidy ploidy, double contamination, double same);

  protected AbstractSomaticCaller getSomaticCaller(final Hypotheses<Description> hypotheses, VariantParams params) {
    return getSomaticCaller(hypotheses, params, 1.0, 1.0);
  }

  protected abstract AbstractSomaticCaller getSomaticCaller(final Hypotheses<Description> hypotheses, VariantParams params, double phi, double psi);

  protected VariantOutputVcfFormatter getFormatter() {
    return getFormatter(null);
  }
  protected VariantOutputVcfFormatter getFormatter(VariantParams params) {
    final VariantOutputVcfFormatter formatter = params == null ?  new VariantOutputVcfFormatter("normal", "cancer") : new VariantOutputVcfFormatter(params, "normal", "cancer");
    formatter.addExtraInfoFields(EnumSet.of(VcfInfoField.LOH, VcfInfoField.NCS));
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
    final double[] qualities = new double[numReads];
    for (int i = 0; i < numReads; i++) {
      qualities[i] = 0.05;
    }
    return doReads(readNt, qualities);
  }

  protected final List<ModelInterface<D>> doReads(final int readNt, double... qualities) {
    return getIncrementer(Ploidy.HAPLOID, 0.0, 0.99).doReads(readNt, qualities).freeze();
  }

  protected List<ModelInterface<Description>> doNormalReads(final int numReads, final int readNt) {
    return getNormalIncremeter(Ploidy.HAPLOID).doReads(numReads, readNt).freeze();
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
      EXPECT_ALL_SAME
    );
  }

  protected static final String EXPECT_NORMAL_EQ_REF = "chr1\t14\t.\tA\tC\t.\tPASS\tNCS=10.864;DP=6\tGT:DP:RE:AR:GQ:ABP:SBP:RPB:PUR:RS:AD:SSC:SS\t0:3:0.293:0.000:36:0.00:6.51:0.00:0.00:A,3,0.293:3,0\t1:3:0.293:0.000:11:0.00:6.51:0.00:0.00:C,3,0.293:0,3:1.0:2\n";

  public void testNormalEqualsRef() throws InvalidParamsException, IOException {
    checkCancer(
      doNormalReads(3, DNARangeAT.A),
      doReads(3, DNARangeAT.C),
      EXPECT_NORMAL_EQ_REF
    );
  }

  protected static final String EXPECT_CANCER_EQ_REF = "chr1\t14\t.\tA\t.\t.\tPASS\tDP=6\tGT:DP:RE:AR:GQ:ABP:SBP:RPB:PUR:RS:AD:SS\t0:3:0.293:0.000:14:26.06:.:0.00:0.00:C,3,0.293:0\t0:3:0.293:0.000:25:0.00:6.51:0.00:0.00:A,3,0.293:3:0\n";

  public void testCancerEqualsRef() throws InvalidParamsException, IOException {
    checkCancer(
      doNormalReads(3, DNARangeAT.C),
      doReads(3, DNARangeAT.A),
      EXPECT_CANCER_EQ_REF
    );
  }

  protected static final String EXPECT_CANCER_EQ_NORMAL = "chr1\t14\t.\tA\tC\t.\tPASS\tDP=6\tGT:DP:RE:AR:GQ:ABP:SBP:RPB:PUR:RS:AD:SS\t1:3:0.293:0.000:55:0.00:6.51:0.00:0.00:C,3,0.293:0,3\t1:3:0.293:0.000:65:0.00:6.51:0.00:0.00:C,3,0.293:0,3:1\n";

  public void testCancerEqualsNormal() throws InvalidParamsException, IOException {
    checkCancer(
      doNormalReads(3, DNARangeAT.C),
      doReads(3, DNARangeAT.C),
      EXPECT_CANCER_EQ_NORMAL
    );
  }

  protected static final String EXPECT_ALL_DIFFERENT = "chr1\t14\t.\tA\tC,G\t.\tPASS\tNCS=8.224;DP=6\tGT:DP:RE:AR:GQ:ABP:SBP:RPB:PUR:RS:AD:SSC:SS\t1:3:0.293:0.000:11:0.00:6.51:0.00:0.00:C,3,0.293:0,3,0\t2:3:0.293:0.000:11:0.00:6.51:0.00:0.00:G,3,0.293:0,0,3:0.7:2\n";

  public void testAllDifferent() throws InvalidParamsException, IOException {
    checkCancer(
      doNormalReads(3, DNARangeAT.C),
      doReads(3, DNARangeAT.G),
      EXPECT_ALL_DIFFERENT
    );
  }

  private void checkCancer(List<ModelInterface<Description>> normal, List<ModelInterface<D>> cancer, String expect) throws InvalidParamsException, IOException {
    checkCancer(normal, cancer, expect, getDefaultParams());
  }

  VariantParams getDefaultParams() {
    final GenomeRelationships genomeRelationships = getGenomeRelationships();
    return new VariantParamsBuilder().callLevel(VariantOutputLevel.ALL).somaticParams(new SomaticParamsBuilder().somaticRate(0.001).create()).genomeRelationships(genomeRelationships).create();
  }

  private void checkCancer(List<ModelInterface<Description>> normal, List<ModelInterface<D>> cancer, String expect, VariantParams params) throws InvalidParamsException, IOException {
    final Variant v = getVariant(normal, cancer, params);
    final VariantOutputVcfFormatter formatter = getFormatter(params);
    assertEquals(expect, formatter.formatCall(v));
  }

  Variant getVariant(List<ModelInterface<Description>> normal, List<ModelInterface<D>> cancer, VariantParams params) {
    return getVariant(normal, cancer, params, 1.0, 1.0);
  }
  Variant getVariant(List<ModelInterface<Description>> normal, List<ModelInterface<D>> cancer, VariantParams params, double phi, double psi) {
    final int refNt = DNARange.A;
    final Hypotheses<Description> hypotheses = normal.get(0).hypotheses();
    final AbstractSomaticCaller ccs = getSomaticCaller(hypotheses, params, phi, psi);
    ccs.integrity();
    //System.out.println(new Posterior(ccs.mQ, A, normal.categories(), cancer.categories()));
    //    final MultivarianceCall out = new MultivarianceCall("chr1", 13, 14, VarianceCallType.SNP);
    //    out.setName("foo");
    final byte[] ref = new byte[14];
    ref[12] = DNARange.G;
    ref[13] = refNt;
    //debugCalls(outParams, 13, normal[refCode], cancer[refCode]);

    final List<ModelInterface<?>> models = new ArrayList<>();
    models.add(normal.get(0));
    models.add(cancer.get(0));
//    for (ModelInterface<?> model : models) {
//      model.freeze();
//    }
    final HypothesesPrior<Description> hypotheses1 = (HypothesesPrior<Description>) normal.get(0).hypotheses();
    final HypothesesPrior<Description> haploid;
    final HypothesesPrior<Description> diploid;
    if (hypotheses1.haploid()) {
      haploid = hypotheses1;
      diploid = null;
    } else {
      haploid = null;
      diploid = hypotheses1;
    }
    final HaploidDiploidHypotheses<HypothesesPrior<Description>> hyp = new HaploidDiploidHypotheses<>(HypothesesNone.SINGLETON,
      haploid,
        diploid);
    return ccs.makeCall("chr1", 13, 14, ref, models, hyp);
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

  protected static final String EXPECT_VAF = "chr1\t14\t.\tA\tG\t.\tPASS\tDP=3\tGT:VA:DP:RE:AR:GQ:ABP:SBP:RPB:PUR:RS:ADE:AD:SS:VAF\t1:.:0:0.000:.:10:.:.:.:.:.:0.0,0.0:0,0\t1:1:3:0.293:0.000:21:0.00:6.51:0.00:0.00:G,3,0.293:0.0,2.7:0,3:1:1.000\n";
  public void testVariantAllele() throws InvalidParamsException, IOException {

    final VariantParams params = new VariantParamsBuilder()
      .callLevel(VariantOutputLevel.ALL)
      .somaticParams(new SomaticParamsBuilder().somaticRate(0.001).create())
      .minVariantAllelicFraction(0.1)
      .minVariantAllelicDepth(1)
      .genomeRelationships(getGenomeRelationships())
      .create();
    checkCancer(
      doNormalReads(0, DNARangeAT.C),
      doReads(3, DNARangeAT.G),
      EXPECT_VAF,
      params
    );
  }

  GenomeRelationships getGenomeRelationships() {
    final GenomeRelationships genomeRelationships = new GenomeRelationships();
    genomeRelationships.addGenome("normal");
    genomeRelationships.addGenome("cancer");
    return genomeRelationships;
  }

  public void testIncrementNormal() {
    final Variant confident = getVariant(
      getNormalIncremeter(Ploidy.DIPLOID, 0.9)
        .doReads(50, DNARangeAT.C)
        .freeze(),
      getIncrementer(Ploidy.DIPLOID, 0.0, 0.99).doReads(50, DNARangeAT.C)
        .doReads(50, DNARangeAT.G)
        .freeze(),
      getDefaultParams());
    Double lastCancerScore = confident.getNormalCancerScore();
    int i;
    for (i = 1; i < 50; i++) {
      final Variant suspect = getVariant(
        getNormalIncremeter(Ploidy.DIPLOID, 0.9)
          .doReads(50, DNARangeAT.C, 0.0001)
          .doReads(i, DNARangeAT.G, 0.0001)
          .freeze(),
        getIncrementer(Ploidy.DIPLOID, 0.0, 0.99).doReads(50, DNARangeAT.C, 0.0001)
          .doReads(50, DNARangeAT.G, 0.0001)
          .freeze(),
        getDefaultParams());

      final Double suspectCancerScore = suspect.getNormalCancerScore();
      if (suspectCancerScore == null) {
        break;
      }
      assertTrue(suspectCancerScore + " < " + lastCancerScore, suspectCancerScore < lastCancerScore);
      lastCancerScore = suspectCancerScore;
    }
  }
  int firstNonCancerCall(double phi, double psi) {
    final Variant confident = getVariant(
      getNormalIncremeter(Ploidy.DIPLOID, 0.9)
        .doReads(50, DNARangeAT.C)
        .freeze(),
      getIncrementer(Ploidy.DIPLOID, 0.0, 0.99).doReads(50, DNARangeAT.C)
        .doReads(50, DNARangeAT.G)
        .freeze(),
      getDefaultParams());
    Double lastCancerScore = confident.getNormalCancerScore();
    int i;
    for (i = 1; i < 50; i++) {
      final Variant suspect = getVariant(
        getNormalIncremeter(Ploidy.DIPLOID, 0.9)
          .doReads(50, DNARangeAT.C, 0.0001)
          .doReads(i, DNARangeAT.G, 0.0001)
          .freeze(),
        getIncrementer(Ploidy.DIPLOID, 0.0, 0.99).doReads(50, DNARangeAT.C, 0.0001)
          .doReads(50, DNARangeAT.G, 0.0001)
          .freeze(),
        getDefaultParams(), phi, psi);

      if (suspect.getSample(1).isDeNovo() != VariantSample.DeNovoStatus.IS_DE_NOVO) {
        break;
      }
      final Double suspectCancerScore = suspect.getNormalCancerScore();
      assertTrue(suspectCancerScore + " < " + lastCancerScore, suspectCancerScore < lastCancerScore);
      lastCancerScore = suspectCancerScore;
    }
    return i;
  }
  public void testContraryEvidenceCausesQuickerNonCancerCalls() {
    final int baseline = firstNonCancerCall(1.0, 1.0);
    final int contraryEvidence = firstNonCancerCall(0.01, 0.01);
    assertTrue(contraryEvidence + " < " + baseline, contraryEvidence < baseline);

  }

  public void testLessEvidenceCancer() {
    final Variant confident = getVariant(
      getNormalIncremeter(Ploidy.DIPLOID, 0.9)
        .doReads(20, DNARangeAT.A)
        .freeze(),
      getIncrementer(Ploidy.DIPLOID, 0.0, 0.99)
        .doReads(10, DNARangeAT.A)
        .doReads(10, DNARangeAT.C)
        .freeze(),
      getDefaultParams());
    final Variant suspect = getVariant(
      getNormalIncremeter(Ploidy.DIPLOID, 0.9)
        .doReads(20, DNARangeAT.A)
        .freeze(),
      getIncrementer(Ploidy.DIPLOID, 0.0, 0.99)
        .doReads(10, DNARangeAT.A)
        .doReads(9, DNARangeAT.C)
        .freeze(),
      getDefaultParams());
    final double confidentCancerScore = confident.getNormalCancerScore();
    final double suspectCancerScore = suspect.getNormalCancerScore();
    assertTrue(suspectCancerScore + " < " + confidentCancerScore, suspectCancerScore < confidentCancerScore);
  }

  public void testMoreEvidenceCancer() {
    final Variant confident = getVariant(
      getNormalIncremeter(Ploidy.DIPLOID, 0.9)
        .doReads(20, DNARangeAT.A)
        .freeze(),
      getIncrementer(Ploidy.DIPLOID, 0.0, 0.99)
        .doReads(10, DNARangeAT.A)
        .doReads(10, DNARangeAT.C)
        .freeze(),
      getDefaultParams());
    final Variant improved = getVariant(
      getNormalIncremeter(Ploidy.DIPLOID, 0.9)
        .doReads(20, DNARangeAT.A)
        .freeze(),
      getIncrementer(Ploidy.DIPLOID, 0.0, 0.99)
        .doReads(10, DNARangeAT.A)
        .doReads(11, DNARangeAT.C)
        .freeze(),
      getDefaultParams());
    final double confidentCancerScore = confident.getNormalCancerScore();
    final double improvedScore = improved.getNormalCancerScore();
    assertTrue(improvedScore + " > " + confidentCancerScore, improvedScore > confidentCancerScore);
  }
}
