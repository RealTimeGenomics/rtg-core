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
    for (int i = 0; i < result.length; i++) {
      result[i] = arith.prob2Poss(vals[i]) / MathUtils.LOG_10;
    }
    return result;
  }
  public void assertApproxEquals(double[] a, double[] b, double tolerance) {
    String failMsg = Arrays.toString(a) + " != " + Arrays.toString(b);
    if (a.length != b.length) {
      fail(failMsg);
    }
    for (int i = 0; i < a.length; i++) {
      assertEquals(failMsg, a[i], b[i], tolerance);
    }
  }
  public void test() {
    final VariantSample sample = getVariantSampleHaploid(0.1, 0.3, 0.4, 0.6);
    final List<String> calls = Arrays.asList("A", "C", "T");
    sample.setGenotypeLikelihoods(sample.computeGenotypeLikelihoods(new HashSet<>(calls)));
    final double[] likelihoods = GenotypeLikelihoodUtils.getLikelihoods(sample, calls);

    final double[] expected = convert(0.1, 0.3, 0.6);
    assertApproxEquals(expected, likelihoods, 0.0001);
  }
  public void testDiploid() {
    final VariantSample sample = getVariantSampleDiploid(0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1);
    final List<String> calls = Arrays.asList("A", "T");
    sample.setGenotypeLikelihoods(sample.computeGenotypeLikelihoods(new HashSet<>(calls)));
    final double[] likelihoods = GenotypeLikelihoodUtils.getLikelihoods(sample, calls);

    final double[] expected = convert(1.0 / 3, 1.0 / 3, 1.0 / 3);
    assertApproxEquals(expected, likelihoods, 0.0001);
  }
  public void testDiploid2() {
    final VariantSample sample = getVariantSampleDiploid(0.2, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1);
    final List<String> calls = Arrays.asList("A", "T");
    sample.setGenotypeLikelihoods(sample.computeGenotypeLikelihoods(new HashSet<>(calls)));
    final double[] likelihoods = GenotypeLikelihoodUtils.getLikelihoods(sample, calls);

    final double[] expected = convert(2.0 / 4, 1.0 / 4, 1.0 / 4);
    assertApproxEquals(expected, likelihoods, 0.0001);
  }
  public void testDiploid3() {
    final VariantSample sample = getVariantSampleDiploid(0.2, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.3);
    final List<String> calls = Arrays.asList("A", "T");
    sample.setGenotypeLikelihoods(sample.computeGenotypeLikelihoods(new HashSet<>(calls)));
    final double[] likelihoods = GenotypeLikelihoodUtils.getLikelihoods(sample, calls);

    final double[] expected = convert(2.0 / 6, 3.0 / 6, 1.0 / 6);
    assertApproxEquals(expected, likelihoods, 0.0001);
  }
}
