/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
