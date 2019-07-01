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

package com.rtg.variant.bayes.snp;

import com.rtg.variant.VariantOutputOptions;
import com.rtg.variant.VariantParams;
import com.rtg.variant.VariantSample;
import com.rtg.variant.bayes.EvidenceInterface;
import com.rtg.variant.bayes.Hypotheses;
import com.rtg.variant.bayes.Statistics;
import com.rtg.variant.util.VariantUtils;
import com.rtg.variant.util.arithmetic.SimplePossibility;

import junit.framework.TestCase;

/**
 */
public class StatisticsSnpTest extends TestCase {

  static final class MockHypotheses extends Hypotheses<DescriptionCommon> {
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

  private static EvidenceInterface di(final int read, final int score, double r, boolean unmapped) {
    return new EvidenceQ(DescriptionSnp.SINGLETON, read, 0, 0, r, VariantUtils.phredToProb(score), true, false, true, false, unmapped);
  }

  //empty state
  public void test0() {
    final Hypotheses<DescriptionCommon> hy = new MockHypotheses(0, false);
    final VariantOutputOptions params = VariantParams.builder().create();
    final Statistics<?> ss = new StatisticsSnp(hy.description());
    assertEquals(0, ss.nonRefCount());
    assertEquals(0, ss.ambiguousCount());
    assertEquals(null, ss.ambiguityRatio());
    assertFalse(ss.overCoverage(params, ""));
    assertFalse(ss.ambiguous(params));

    final VariantSample vs = new VariantSample(null);
    ss.addCountsToSample(vs, null, params);
    assertEquals("", vs.getStatisticsString().replace("\t", " "));
    assertEquals(Integer.valueOf(0), vs.getCoverage());
    assertEquals(0.0, vs.getCorrection());
  }

  //nonRef and also reset followed by increment
  public void test1() {
    final Hypotheses<DescriptionCommon> hy = new MockHypotheses(0, false);
    final Statistics<?> ss = new StatisticsSnp(hy.description());

    final EvidenceInterface di = di(1, 10, 0.2, false);
    ss.increment(di, hy.reference());
    assertEquals(0, ss.placedUnmappedCount());

    final VariantOutputOptions params = VariantParams.builder().create();
    final VariantSample vs = new VariantSample(null);
    ss.addCountsToSample(vs, null, params);
    assertEquals(" C 1 0.280", vs.getStatisticsString().replace("\t", " "));
    assertEquals(0.280, vs.getCorrection());
  }

  //nothing special
  public void test2() {
    final Hypotheses<DescriptionCommon> hy = new MockHypotheses(0, false);
    final Statistics<?> ss = new StatisticsSnp(hy.description());

    final EvidenceInterface di = di(0, 10, 0.2, false);
    ss.increment(di, hy.reference());

    final VariantOutputOptions params = VariantParams.builder().create();
    final VariantSample vs = new VariantSample(null);
    ss.addCountsToSample(vs, null, params);
    assertEquals(" A 1 0.280", vs.getStatisticsString().replace("\t", " "));
    assertEquals(0.280, vs.getCorrection());
  }

  public void testAmbiguous() {
    final Hypotheses<DescriptionCommon> hy = new MockHypotheses(0, false);
    final Statistics<?> ss = new StatisticsSnp(hy.description());

    final EvidenceInterface di = di(0, 0, 0.2, false);
    ss.increment(di, hy.reference());

    final VariantOutputOptions params = VariantParams.builder().create();
    final VariantSample vs = new VariantSample(null);
    ss.addCountsToSample(vs, null, params);
    assertEquals(" A 1 1.000", vs.getStatisticsString().replace("\t", " "));

    assertEquals(1.0, vs.getCorrection());
  }

  public void testPlacedUnmapped() {
    final Hypotheses<DescriptionCommon> hy = new MockHypotheses(0, false);
    final Statistics<?> ss = new StatisticsSnp(hy.description());
    final EvidenceInterface di = di(0, 0, 0.2, true);
    ss.increment(di, hy.reference());
    assertEquals(1, ss.placedUnmappedCount());
  }

}
