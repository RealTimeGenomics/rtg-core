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

import com.rtg.reference.Ploidy;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.variant.StaticThreshold;
import com.rtg.variant.Variant;
import com.rtg.variant.Variant.VariantFilter;
import com.rtg.variant.VariantLocus;
import com.rtg.variant.VariantOutputOptions;
import com.rtg.variant.VariantParams;
import com.rtg.variant.VariantSample;
import com.rtg.variant.bayes.snp.DescriptionCommon;
import com.rtg.variant.bayes.snp.DescriptionNone;
import com.rtg.variant.bayes.snp.DescriptionSnp;
import com.rtg.variant.bayes.snp.EvidenceQ;
import com.rtg.variant.util.VariantUtils;
import com.rtg.variant.util.arithmetic.SimplePossibility;

import junit.framework.TestCase;

/**
 */
public abstract class AbstractStatisticsTest extends TestCase {

  private static final class MockHypotheses extends Hypotheses<DescriptionCommon> {
    private final int mRef;

    MockHypotheses(final int ref, final boolean haploid) {
      super(DescriptionSnp.SINGLETON, SimplePossibility.SINGLETON, haploid);
      mRef = ref;
    }
    @Override
    public int reference() {
      return mRef;
    }
  }

  protected abstract Statistics<?> getStatistics(Description d);

  private void incrementStatistics(Statistics<?> stats, Hypotheses<DescriptionCommon> hyp, EvidenceInterface distribution) {
    stats.increment(distribution, hyp.reference());
  }

  static EvidenceInterface di(final int read, final int score, double r, boolean forward, boolean pairedRead, boolean mated) {
    return new EvidenceQ(DescriptionSnp.SINGLETON, read, 0, 0, r, VariantUtils.phredToProb(score), forward, pairedRead, true, mated, false);
  }

  public void testAmbiguity1() {
    checkAmbiguity(0.9999, true);
    checkAmbiguity(1.0001, false);
    checkAmbiguity(null, false);
  }

  private void checkAmbiguity(final Double maxAmb, final boolean amb) {
    final VariantOutputOptions params = VariantParams.builder().maxAmbiguity(maxAmb).create();
    final Hypotheses<DescriptionCommon> hy = new MockHypotheses(0, false);
    final Statistics<?> ss = getStatistics(hy.description());

    final EvidenceInterface di = di(0, 0, 0.51, true, false, false);
    incrementStatistics(ss, hy, di);
    incrementStatistics(ss, hy, di);
    assertEquals(2, ss.ambiguousCount());
    assertEquals(1.0, ss.ambiguityRatio());
    assertEquals(amb, ss.ambiguous(params));

    final Variant v = new Variant(new VariantLocus("foo", 1, 2, "c", (char) -1));
    ss.addFiltersToVariant(v, null, params);
    final VariantSample vs = new VariantSample(null);
    ss.addCountsToSample(vs, null, params);
    assertEquals(amb, v.isFiltered(VariantFilter.AMBIGUITY));
    assertEquals(1.0, vs.getAmbiguityRatio());
  }

  public void testAlleleBalance() {
    final VariantOutputOptions params = VariantParams.builder().create();
    final Hypotheses<DescriptionCommon> hy = new MockHypotheses(0, false);
    final Statistics<?> ss = getStatistics(hy.description());

    incrementStatistics(ss, hy, di(0, 10, 0.0, true, false, false));
    incrementStatistics(ss, hy, di(1, 10, 0.0, true, false, false));

    final VariantSample vs = new VariantSample(null);
    ss.addCountsToSample(vs, null, params);
  }

  public void testCoverage() {
    checkCoverage(true, 1);
    checkCoverage(false, 2);
  }

  private void checkCoverage(final boolean coverage, final int threshold) {
    final VariantOutputOptions params = VariantParams.builder().maxCoverageFilter(new StaticThreshold(threshold)).create();
    final Hypotheses<DescriptionCommon> hy = new MockHypotheses(0, false);
    final Statistics<?> ss = getStatistics(hy.description());

    incrementStatistics(ss, hy, di(0, 10, 0.0, true, false, false));
    incrementStatistics(ss, hy, di(1, 10, 0.0, true, false, false));
    assertEquals(coverage, ss.overCoverage(params, ""));

    final Variant v = new Variant(new VariantLocus("foo", 1, 2, "c", (char) -1));
    ss.addFiltersToVariant(v, null, params);
    assertEquals(coverage, v.isFiltered(VariantFilter.COVERAGE));
  }

  //empty state
  public void test0() {
    final VariantOutputOptions params = VariantParams.builder().create();
    final Statistics<?> ss = getStatistics(DescriptionNone.SINGLETON);
    assertEquals(0, ss.coverage());
    assertEquals(0, ss.nonRefCount());
    assertEquals(0, ss.ambiguousCount());
    assertEquals(null, ss.ambiguityRatio());
    assertFalse(ss.overCoverage(params, ""));
    assertFalse(ss.ambiguous(params));

    final VariantSample sample = new VariantSample(null);
    ss.addCountsToSample(sample, null, params);
    assertNull(sample.getStatisticsString());
    assertEquals(Integer.valueOf(0), sample.getCoverage());
    assertEquals((double) 0, sample.getCorrection());
  }

  //nonRef and also reset followed by increment
  public void test1() {
    final Hypotheses<DescriptionCommon> hy = new MockHypotheses(0, false);
    final Statistics<?> ss = getStatistics(hy.description());

    final EvidenceInterface di = di(1, 10, 0.2, true, false, false);
    incrementStatistics(ss, hy, di);
    assertEquals(1, ss.coverage());
    assertEquals(1, ss.nonRefCount());
    assertEquals(0, ss.ambiguousCount());
    final VariantOutputOptions params = VariantParams.builder().create();
    final VariantSample sample = new VariantSample(null);
    ss.addCountsToSample(sample, null, params);
    assertNull(sample.getStatisticsString());
    assertEquals(Integer.valueOf(1), sample.getCoverage());
    assertEquals(0.28, sample.getCorrection());
  }

  //nothing special
  public void test2() {
    final Hypotheses<DescriptionCommon> hy = new MockHypotheses(0, false);
    final Statistics<?> ss = getStatistics(hy.description());

    final EvidenceInterface di = di(0, 10, 0.2, true, false, false);
    incrementStatistics(ss, hy, di);
    assertEquals(1, ss.coverage());
    assertEquals(0, ss.nonRefCount());
    assertEquals(0, ss.ambiguousCount());

    final VariantOutputOptions params = VariantParams.builder().create();
    final VariantSample sample = new VariantSample(null);
    ss.addCountsToSample(sample, null, params);
    assertNull(sample.getStatisticsString());
    assertEquals(Integer.valueOf(1), sample.getCoverage());
    assertEquals(0.28, sample.getCorrection());
  }

  //ambiguous
  public void test3() {
    final Hypotheses<DescriptionCommon> hy = new MockHypotheses(0, false);
    final Statistics<?> ss = getStatistics(hy.description());

    final EvidenceInterface di = di(0, 0, 0.51, true, false, false);
    incrementStatistics(ss, hy, di);
    assertEquals(1, ss.coverage());
    assertEquals(0, ss.nonRefCount());
    assertEquals(1, ss.ambiguousCount());

    final VariantOutputOptions params = VariantParams.builder().create();
    final VariantSample sample = new VariantSample(null);
    ss.addCountsToSample(sample, null, params);
    assertNull(sample.getStatisticsString());
    assertEquals(Integer.valueOf(1), sample.getCoverage());
    assertEquals(1.0, sample.getCorrection());
  }

  public void testStrandBias() {
    final Hypotheses<DescriptionCommon> hy = new MockHypotheses(0, false);
    final Statistics<?> ss = getStatistics(hy.description());

    final EvidenceInterface di = di(0, 0, 0.51, true, false, false);
    incrementStatistics(ss, hy, di);
    incrementStatistics(ss, hy, di);
    final EvidenceInterface di2 = di(0, 0, 0.51, false, false, false);
    incrementStatistics(ss, hy, di2);
  }
  public void testUnmated() {
    final Hypotheses<DescriptionCommon> hy = new MockHypotheses(1, false);
    final Statistics<?> ss = getStatistics(hy.description());

    final EvidenceInterface di = di(0, 0, 0.51, false, true, true);
    incrementStatistics(ss, hy, di);
    incrementStatistics(ss, hy, di);
    final EvidenceInterface di2 = di(0, 0, 0.51, false, true, false);
    incrementStatistics(ss, hy, di2);
    final EvidenceInterface di3 = di(1, 0, 0.51, false, true, false);
    incrementStatistics(ss, hy, di3);
    final VariantSample v = new VariantSample(Ploidy.HAPLOID, "A", false, new MockGenotypeMeasure(1.0), VariantSample.DeNovoStatus.UNSPECIFIED, null);
    ss.addCountsToSample(v, new MockModel<>(hy, ss, new double[hy.size()]), VariantParams.builder().create());
    assertEquals(0.7238, v.getHoeffdingUnmatedBiasAllele1(), 0.0001);
    assertNull(v.getHoeffdingStrandBiasAllele2());
  }
  public void testUnmatedDiploid() {
    final Hypotheses<DescriptionCommon> hy = new MockHypotheses(1, false);
    final Statistics<?> ss = getStatistics(hy.description());

    final EvidenceInterface di = di(0, 0, 0.51, false, true, true);
    incrementStatistics(ss, hy, di);
    incrementStatistics(ss, hy, di);
    final EvidenceInterface di2 = di(0, 0, 0.51, false, true, false);
    incrementStatistics(ss, hy, di2);
    final EvidenceInterface di3 = di(1, 0, 0.51, false, true, false);
    incrementStatistics(ss, hy, di3);
    incrementStatistics(ss, hy, di3);
    incrementStatistics(ss, hy, di3);
    incrementStatistics(ss, hy, di3);
    incrementStatistics(ss, hy, di3);
    final EvidenceInterface di4 = di(1, 0, 0.51, false, true, true);
    incrementStatistics(ss, hy, di4);
    final VariantSample v = new VariantSample(Ploidy.DIPLOID, "A" + VariantUtils.COLON + "C", false, new MockGenotypeMeasure(1.0), VariantSample.DeNovoStatus.UNSPECIFIED, null);
    ss.addCountsToSample(v, new MockModel<>(hy, ss, new double[hy.size()]), VariantParams.builder().create());
    assertEquals(2.8952, v.getHoeffdingUnmatedBiasAllele1(), 0.0001);
    assertEquals(1.4476, v.getHoeffdingUnmatedBiasAllele2(), 0.0001);
  }

  public void testAmbiguityRatio() {
    Diagnostic.setLogStream();
    final Hypotheses<DescriptionCommon> hy = new MockHypotheses(1, false);
    final Statistics<?> stat = getStatistics(hy.description());
    for (int i = 0; i < 3; ++i) {
      //mapq = 0
      stat.increment(di(0, 0, VariantUtils.phredToProb(0), true, true, true), 0);
    }
    for (int i = 0; i < 10; ++i) {
      //mapq 20
      stat.increment(di(2, 0, VariantUtils.phredToProb(20), true, true, true), 0);
    }
    final VariantOutputOptions params = VariantParams.builder().create();
    final VariantSample vs = new VariantSample(null, null, true, null, VariantSample.DeNovoStatus.UNSPECIFIED, null);
    stat.addCountsToSample(vs, null, params);
    assertEquals(3 / 13.0, stat.ambiguityRatio());
  }
}
