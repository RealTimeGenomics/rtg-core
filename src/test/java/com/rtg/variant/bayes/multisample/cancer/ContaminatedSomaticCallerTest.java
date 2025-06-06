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
package com.rtg.variant.bayes.multisample.cancer;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import com.rtg.mode.DNA;
import com.rtg.reference.Ploidy;
import com.rtg.relation.GenomeRelationships;
import com.rtg.relation.Relationship;
import com.rtg.variant.SomaticParamsBuilder;
import com.rtg.variant.Variant;
import com.rtg.variant.VariantOutputLevel;
import com.rtg.variant.VariantParams;
import com.rtg.variant.VariantParamsBuilder;
import com.rtg.variant.VariantSample;
import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.Hypotheses;
import com.rtg.variant.bayes.ModelInterface;
import com.rtg.variant.bayes.NoAlleleBalance;
import com.rtg.variant.bayes.snp.StatisticsSnp;
import com.rtg.variant.util.arithmetic.LogPossibility;

/**
 */
public class ContaminatedSomaticCallerTest extends AbstractSomaticCallerTest<Description> {

  @Override
  protected AbstractSomaticCaller getSomaticCaller(Hypotheses<Description> hypotheses, VariantParams params, double phi, double psi) {
    return new ContaminatedSomaticCaller(new DefaultSomaticPriorsFactory<>(hypotheses, 0.0), new DefaultSomaticPriorsFactory<>(hypotheses, 0.0), params, phi, psi, 0.5);
  }

  @Override
  protected List<ModelInterface<Description>> getModel(Ploidy ploidy, double contamination, double same) {
    final List<ModelInterface<Description>> models = new ArrayList<>();
    for (int ref = 0; ref < 4; ++ref) {
      final Hypotheses<Description> hyps = simpleHyps(0.99, ref, ploidy);
      final HypothesesCancer<Hypotheses<Description>> hypc = new HypothesesCancer<>(hyps, LogPossibility.SINGLETON);
      models.add(new ModelCancerContamination<>(hypc, contamination, new StatisticsSnp(hyps.description()), new NoAlleleBalance()));
    }
    return models;
  }

  private Variant getVariant(double contamination, int[] readCounts, int[] normalReadCounts) {
    final GenomeRelationships genomeRelationships = new GenomeRelationships();
    final Relationship relationship = new Relationship("normal", "cancer", Relationship.RelationshipType.ORIGINAL_DERIVED);
    relationship.setProperty(Relationship.CONTAMINATION, String.valueOf(contamination));
    relationship.setProperty(Relationship.REVERSE_CONTAMINATION, String.valueOf(0.001));
    genomeRelationships.addRelationship(relationship);
    final VariantParams variantParams = new VariantParamsBuilder().callLevel(VariantOutputLevel.ALL).somaticParams(new SomaticParamsBuilder().lohPrior(0.1).includeGainOfReference(true).somaticRate(0.000165).create()).genomeRelationships(genomeRelationships).create();
    final ModelIncrementer<Description> incrementer = getIncrementer(Ploidy.DIPLOID, contamination, 0.999);
    for (int i = 0; i < 4; ++i) {
      for (int j = 0; j < readCounts[i]; ++j) {
        incrementer.doRead(DNA.values()[i + 1].ordinal() - 1, 0.005, 0.005);
      }
    }
    final ModelIncrementer<Description> normalIncremeter = getNormalIncremeter(Ploidy.DIPLOID, 0.999);
    for (int i = 0; i < 4; ++i) {
      for (int j = 0; j < normalReadCounts[i]; ++j) {
        normalIncremeter.doRead(DNA.values()[i + 1].ordinal() - 1, 0.05, 0.05);
      }
    }
    return getVariant(
      normalIncremeter
        .freeze(),
      incrementer.freeze(),
      variantParams);

  }

  private void checkTumorOnly(String expectedNormal, String expectedCancer, double contamination, int[] readCounts) {
    checkTumorOnly(expectedNormal, expectedCancer, contamination, readCounts, new int[] {0, 0, 0, 0});
  }

  private void checkTumorOnly(String expectedNormal, String expectedCancer, double contamination, int[] readCounts, int[] normalReadCounts) {
    final Variant call = getVariant(contamination, readCounts, normalReadCounts);
    final String normalName = call.getSample(0).getName();
    final String tumorName = call.getSample(1).getName();
    final String failure = "contamination: " + contamination + ", normal:" + normalName + ", tumor:" + tumorName + ", reads: " + Arrays.toString(readCounts);
//    System.err.println(Arrays.toString(readCounts) + " " + normalName + "_" + tumorName + " " + call.getNormalCancerScore());
    assertEquals(failure, expectedNormal + ", " + expectedCancer, normalName + ", " + tumorName);
  }

  private void assertRange(List<Variant> variants, String normal, String cancer, int start, int end) {
    for (int i = start; i < end; ++i) {
      final Variant v = variants.get(i);
      final String normalName = v.getSample(0).getName();
      final String tumorName = v.getSample(1).getName();
      final String failure = "normal:" + normalName + ", tumor:" + tumorName + ", reads: " + i;
      assertTrue(failure, normalName.equals(normal) && tumorName.equals(cancer));
    }
  }

  private void rangeTest(List<String> normal, List<String> tumor, List<Integer> breakpoints, List<Variant> variants) {
    int start = 0;
    for (int i = 0; i < normal.size(); ++i) {
      final int breakpoint = breakpoints.get(i);
      final String n = normal.get(i);
      final String t = tumor.get(i);
      assertRange(variants, n, t, start, breakpoint);
      start = breakpoint;
    }
  }

  public void testTumorOnly20() {
    final List<Variant> variants = variantRange(0.2);
    rangeTest(
      Arrays.asList("A:A", "A:C", "A:A", "A:C", "A:A", "A:C", "C:C"),
      Arrays.asList("A:A", "A:A", "A:C", "A:C", "C:C", "C:C", "C:C"),
      Arrays.asList(7, 18, 43, 75, 84, 96, 100),
      variants
    );
  }

  public void testTumorOnly30EdgeCase() {
    final Variant variant = getVariant(0.3, new int[]{24, 76, 0, 0}, new int[]{0, 0, 0, 0});
    assertEquals("C:C", variant.getSample(0).getName());
    assertEquals("C:C", variant.getSample(1).getName());
    assertEquals(VariantSample.DeNovoStatus.NOT_DE_NOVO, variant.getSample(1).isDeNovo());
  }

  public void testTumorOnly30() {
    final List<Variant> variants = variantRange(0.3);
    rangeTest(
      Arrays.asList("A:A", "A:C", "A:A", "A:C", "C:C", "C:C", "A:C", "C:C"),
      Arrays.asList("A:A", "A:A", "A:C", "A:C", "A:C", "C:C", "C:C", "C:C"),
      Arrays.asList(8, 17, 41, 73, 76, 77, 95, 100),
      variants
    );
  }

  public void testTumorOnly50() {
    checkTumorOnly("A:A", "A:A", 0.5, new int[]{93, 7, 0, 0});
    checkTumorOnly("A:A", "A:C", 0.5, new int[]{75, 25, 0, 0});
    checkTumorOnly("A:C", "A:C", 0.5, new int[]{60, 40, 0, 0});
  }

  private List<Variant> variantRange(double contamination) {
    return variantRange(contamination, 100);
  }
  private List<Variant> variantRange(double contamination, int numReads) {
    final List<Variant> variants = new ArrayList<>();
    for (int i = 0; i < numReads; ++i) {
      variants.add(getVariant(contamination, new int[] {numReads - i, i, 0, 0}, new int[] {0, 0, 0, 0}));
    }
    return variants;
  }
  public void testTumorOnly60() {
    checkTumorOnly("A:A", "A:C", 0.6, new int[] {80, 20, 0, 0});
  }

  public void testTumorOnly70() {
    checkTumorOnly("A:A", "A:C", 0.7, new int[] {85, 15, 0, 0});
  }

  public void testTumorOnly80() {
    //Too close to read error
    checkTumorOnly("A:A", "A:C", 0.8, new int[] {90, 10, 0, 0});
    // Adding more reads reveals the tumor
    checkTumorOnly("A:A", "A:C", 0.8, new int[] {900, 100, 0, 0});
  }
}
