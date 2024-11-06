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
package com.rtg.variant.format;

import java.util.Arrays;
import java.util.HashSet;
import java.util.List;

import com.rtg.reference.Ploidy;
import com.rtg.util.MathUtils;
import com.rtg.variant.GenomePriorParams;
import com.rtg.variant.VariantSample;
import com.rtg.variant.bayes.ArrayGenotypeMeasure;
import com.rtg.variant.bayes.snp.HypothesesSnp;
import com.rtg.variant.util.arithmetic.LogPossibility;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;
import com.rtg.variant.util.arithmetic.SimplePossibility;

import junit.framework.TestCase;

/**
 */
public class GenotypeLikelihoodUtilsTest extends TestCase {
  private VariantSample getVariantSampleDiploid(double... likely) {
    final HypothesesSnp hypotheses = new HypothesesSnp(SimplePossibility.SINGLETON, GenomePriorParams.builder().create(), false, 0);
    final ArrayGenotypeMeasure measure = new ArrayGenotypeMeasure(SimplePossibility.SINGLETON, likely, hypotheses);
    return new VariantSample(Ploidy.DIPLOID, "A:C", false, measure, VariantSample.DeNovoStatus.NOT_DE_NOVO, 0.0);
  }
  private VariantSample getVariantSampleHaploid(double... likely) {
    final HypothesesSnp hypotheses = new HypothesesSnp(SimplePossibility.SINGLETON, GenomePriorParams.builder().create(), true, 0);
    final ArrayGenotypeMeasure measure = new ArrayGenotypeMeasure(SimplePossibility.SINGLETON, likely, hypotheses);
    return new VariantSample(Ploidy.HAPLOID, "C", false, measure, VariantSample.DeNovoStatus.NOT_DE_NOVO, 0.0);
  }
  public double[] convert(double... vals) {
    final PossibilityArithmetic arith = LogPossibility.SINGLETON;
    final double[] result = new double[vals.length];
    for (int i = 0; i < result.length; ++i) {
      result[i] = arith.prob2Poss(vals[i]) / MathUtils.LOG_10;
    }
    return result;
  }
  public void assertApproxEquals(double[] a, double[] b, double tolerance) {
    final String failMsg = Arrays.toString(a) + " != " + Arrays.toString(b);
    if (a.length != b.length) {
      fail(failMsg);
    }
    for (int i = 0; i < a.length; ++i) {
      assertEquals(failMsg, a[i], b[i], tolerance);
    }
  }

  public void test() {
    final VariantSample sample = getVariantSampleHaploid(0.1, 0.3, 0.4, 0.6);
    final List<String> calls = Arrays.asList("A", "C", "T");
    final double[] likelihoods = GenotypeLikelihoodUtils.haploidLikelihoods(calls, sample.computeGenotypeLikelihoods(new HashSet<>(calls)));

    final double[] expected = convert(0.1, 0.3, 0.6);
    assertApproxEquals(expected, likelihoods, 0.0001);
  }
  public void testDiploid() {
    final VariantSample sample = getVariantSampleDiploid(0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1);
    final List<String> calls = Arrays.asList("A", "T");
    final double[] likelihoods = GenotypeLikelihoodUtils.diploidLikelihoods(calls, sample.computeGenotypeLikelihoods(new HashSet<>(calls)));

    final double[] expected = convert(1.0 / 3, 1.0 / 3, 1.0 / 3);
    assertApproxEquals(expected, likelihoods, 0.0001);
  }
  public void testDiploid2() {
    final VariantSample sample = getVariantSampleDiploid(0.2, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1);
    final List<String> calls = Arrays.asList("A", "T");
    final double[] likelihoods = GenotypeLikelihoodUtils.diploidLikelihoods(calls, sample.computeGenotypeLikelihoods(new HashSet<>(calls)));

    final double[] expected = convert(2.0 / 4, 1.0 / 4, 1.0 / 4);
    assertApproxEquals(expected, likelihoods, 0.0001);
  }
  public void testDiploid3() {
    final VariantSample sample = getVariantSampleDiploid(0.2, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.3);
    final List<String> calls = Arrays.asList("A", "T");
    final double[] likelihoods = GenotypeLikelihoodUtils.diploidLikelihoods(calls, sample.computeGenotypeLikelihoods(new HashSet<>(calls)));

    final double[] expected = convert(2.0 / 6, 3.0 / 6, 1.0 / 6);
    assertApproxEquals(expected, likelihoods, 0.0001);
  }
}
