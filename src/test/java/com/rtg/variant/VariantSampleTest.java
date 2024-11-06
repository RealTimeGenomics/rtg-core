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

package com.rtg.variant;

import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import com.rtg.reference.Ploidy;
import com.rtg.variant.bayes.ArrayGenotypeMeasure;
import com.rtg.variant.bayes.snp.DescriptionSnp;
import com.rtg.variant.bayes.snp.HypothesesSnp;
import com.rtg.variant.bayes.snp.StatisticsSnp;
import com.rtg.variant.util.arithmetic.SimplePossibility;

import junit.framework.TestCase;

/**
 */
public class VariantSampleTest extends TestCase {

  public void testSetters() {
    final VariantSample vs = new VariantSample(Ploidy.HAPLOID, "A:C", false, null, VariantSample.DeNovoStatus.UNSPECIFIED, null);

    assertEquals(Ploidy.HAPLOID, vs.getPloidy());
    assertEquals("A:C", vs.getName());
    assertFalse(vs.isIdentity());
    assertEquals(VariantSample.DeNovoStatus.UNSPECIFIED, vs.isDeNovo());
    assertEquals(0.0, vs.getPosterior());


    vs.setStatisticsString("ABCD");
    assertEquals("A:C", vs.toString());
    assertEquals("ABCD", vs.getStatisticsString());
  }

  public void testCopy() {
    final VariantSample vs = new VariantSample(Ploidy.HAPLOID, "A:C", false, null, VariantSample.DeNovoStatus.NOT_DE_NOVO, 0.0);

    vs.setAmbiguityRatio(0.123);
    vs.setCoverage(25, 0.023);
    vs.setStatisticsString("this is a long statistic string");

    final VariantSample vs2 = new VariantSample(Ploidy.HAPLOID, ":C", false, null, VariantSample.DeNovoStatus.NOT_DE_NOVO, 0.0);

    assertNull(vs2.getAmbiguityRatio());
    assertNull(vs2.getCoverage());
    assertNull(vs2.getCorrection());
    assertNull(vs2.getNonIdentityPosterior());
    assertNull(vs2.getStatisticsString());
    vs.setStats(new StatisticsSnp(DescriptionSnp.SINGLETON));
    VariantSample.copy(vs, vs2);

    assertEquals(0.123, vs2.getAmbiguityRatio());
    assertEquals((Integer) 25, vs2.getCoverage());
    assertEquals(0.023, vs2.getCorrection());
    assertEquals("this is a long statistic string", vs2.getStatisticsString());

  }

  public void testGenotypeLikelihoodA() {
    final VariantSample vs = getVariantSampleDiploid();
    final Map<Set<String>, Double> justA = vs.computeGenotypeLikelihoods(Collections.singleton("A"));
    assertEquals(1, justA.size());
    assertEquals(Math.log(1.0), justA.get(Collections.singleton("A")));
  }
  Set<String> set(String ... alleles) {
    return new HashSet<>(Arrays.asList(alleles));
  }
  public void testGenotypeLikelihoodAC() {
    final VariantSample vs = getVariantSampleDiploid();
    final Set<String> alleles = new HashSet<>();
    alleles.add("A");
    alleles.add("C");
    final Map<Set<String>, Double> likelihoods = vs.computeGenotypeLikelihoods(alleles);
    assertEquals(3, likelihoods.size());
    assertEquals(Math.log(1.0), likelihoods.get(Collections.singleton("A")));
    assertEquals(Math.log(5.0), likelihoods.get(set("A", "C")));
  }
  public void testGenotypeLikelihoodAGT() {
    final VariantSample vs = getVariantSampleDiploid();
    final Set<String> alleles = new HashSet<>();
    alleles.add("A");
    alleles.add("G");
    alleles.add("T");
    final Map<Set<String>, Double> likelihoods = vs.computeGenotypeLikelihoods(alleles);
    assertEquals(6, likelihoods.size());
    assertEquals(Math.log(1.0), likelihoods.get(Collections.singleton("A")));
    assertNull(likelihoods.get(set("A", "C")));
    assertEquals(Math.log(8.0), likelihoods.get(set("A", "G")));
    assertEquals(Math.log(8.0), likelihoods.get(set("G", "A")));
    assertEquals(Math.log(7.0), likelihoods.get(set("T", "G")));
    assertEquals(Math.log(10.0), likelihoods.get(set("A", "T")));
    assertEquals(Math.log(4.0), likelihoods.get(set("T")));
  }

  private VariantSample getVariantSampleDiploid() {
    final HypothesesSnp hypotheses = new HypothesesSnp(SimplePossibility.SINGLETON, GenomePriorParams.builder().create(), false, 0);
    final ArrayGenotypeMeasure measure = new ArrayGenotypeMeasure(SimplePossibility.SINGLETON, new double[] {1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, hypotheses);
    return new VariantSample(Ploidy.DIPLOID, "A:C", false, measure, VariantSample.DeNovoStatus.NOT_DE_NOVO, 0.0);
  }
  private VariantSample getVariantSampleHaploid() {
    final HypothesesSnp hypotheses = new HypothesesSnp(SimplePossibility.SINGLETON, GenomePriorParams.builder().create(), true, 0);
    final ArrayGenotypeMeasure measure = new ArrayGenotypeMeasure(SimplePossibility.SINGLETON, new double[] {1, 2, 3, 4}, hypotheses);
    return new VariantSample(Ploidy.HAPLOID, "A:C", false, measure, VariantSample.DeNovoStatus.NOT_DE_NOVO, 0.0);
  }
  public void testHaploidLikelihood() {
    final VariantSample vs = getVariantSampleHaploid();
    final Set<String> alleles = new HashSet<>();
    alleles.add("A");
    alleles.add("G");
    alleles.add("T");
    final Map<Set<String>, Double> likelihoods = vs.computeGenotypeLikelihoods(alleles);
    assertEquals(3, likelihoods.size());
    assertEquals(Math.log(1.0), likelihoods.get(set("A")));
    assertEquals(Math.log(3.0), likelihoods.get(set("G")));
    assertEquals(Math.log(4.0), likelihoods.get(set("T")));
    assertNull(likelihoods.get(set("C")));
  }

}
