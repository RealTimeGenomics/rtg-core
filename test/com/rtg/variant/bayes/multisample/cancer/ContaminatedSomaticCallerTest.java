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
    return new ContaminatedSomaticCaller(new SomaticPriorsFactory<>(hypotheses, 0.0), new SomaticPriorsFactory<>(hypotheses, 0.0), params, phi, psi);
  }

  @Override
  protected List<ModelInterface<Description>> getModel(Ploidy ploidy, double contamination, double same) {
    final List<ModelInterface<Description>> models = new ArrayList<>();
    for (int ref = 0; ref < 4; ref++) {
      final Hypotheses<Description> hyps = simpleHyps(0.99, ref, ploidy);
      final HypothesesCancer<Hypotheses<Description>> hypc = new HypothesesCancer<>(hyps, LogPossibility.SINGLETON);
      models.add(new ModelCancerContamination<>(hypc, contamination, new StatisticsSnp(hyps.description()), new NoAlleleBalance()));
    }
    return models;
  }

  private Variant tumorOnlyVariant(double contamination, int[] readCounts) {
    final GenomeRelationships genomeRelationships = new GenomeRelationships();
    final Relationship relationship = new Relationship("normal", "cancer", Relationship.RelationshipType.ORIGINAL_DERIVED);
    relationship.setProperty(Relationship.CONTAMINATION, String.valueOf(contamination));
    relationship.setProperty(Relationship.REVERSE_CONTAMINATION, String.valueOf(0.001));
    genomeRelationships.addRelationship(relationship);
    final VariantParams variantParams = new VariantParamsBuilder().callLevel(VariantOutputLevel.ALL).somaticParams(new SomaticParamsBuilder().somaticRate(0.000165).create()).genomeRelationships(genomeRelationships).create();
    final ModelIncrementer<Description> incrementer = getIncrementer(Ploidy.DIPLOID, contamination, 0.999);
    for (int i = 0; i < 4; i++) {
      incrementer.doReads(readCounts[i], DNA.values()[i + 1].ordinal() - 1);
    }
    return getVariant(
      getNormalIncremeter(Ploidy.DIPLOID, 0.999)
        .freeze(),
        incrementer.freeze(),
      variantParams);
  }


  private void checkTumorOnly(String expectedNormal, String expectedCancer, double contamination, int[] readCounts) {
    final Variant call = tumorOnlyVariant(contamination, readCounts);
    final String normalName = call.getSample(0).getName();
    final String tumorName = call.getSample(1).getName();
    final String failure = "contamination: " + contamination + ", normal:" + normalName + ", tumor:" + tumorName + ", reads: " + Arrays.toString(readCounts);
//    System.err.println(Arrays.toString(readCounts) + " " + normalName + "_" + tumorName + " " + call.getNormalCancerScore());
    assertEquals(failure, expectedNormal + ", " + expectedCancer, normalName + ", " + tumorName);
  }

  public void testTumorOnly20() {
    // This is too close to being heterozygous to detect at perfect
    checkTumorOnly("A:C", "A:C", 0.2, new int[] {60, 40, 0, 0});

    // but we can pick up the lower side of the expected VAF
    checkTumorOnly("A:A", "A:C", 0.2, new int[] {65, 35, 0, 0});
  }


  public void testTumorOnly25() {
    checkTumorOnly("A:A", "A:A", 0.25, new int[] {90, 10, 0, 0});
    checkTumorOnly("A:A", "A:C", 0.25, new int[] {65, 35, 0, 0});
    checkTumorOnly("A:C", "A:C", 0.25, new int[] {55, 45, 0, 0});
  }
  public void testTumorOnly30() {
    checkTumorOnly("A:A", "A:A", 0.3, new int[] {90, 10, 0, 0});
    checkTumorOnly("A:A", "A:C", 0.3, new int[] {65, 35, 0, 0});
    checkTumorOnly("A:C", "A:C", 0.3, new int[] {55, 45, 0, 0});
  }
  public void testTumorOnly50() {
    checkTumorOnly("A:A", "A:A", 0.5, new int[]{90, 10, 0, 0});
    checkTumorOnly("A:A", "A:C", 0.5, new int[]{75, 25, 0, 0});
    checkTumorOnly("A:C", "A:C", 0.5, new int[]{60, 40, 0, 0});
  }

  public void rangeCheck(double contamination, int minTumor, int maxTumor) {
    int i = 0;
    for (; i < minTumor; i++) {
      checkTumorOnly("A:A", "A:A", contamination, new int[] {100 - i, i, 0, 0});
    }
    for (; i < maxTumor; i++) {
      checkTumorOnly("A:A", "A:C", contamination, new int[] {100 - i, i, 0, 0});
    }
    for (; i < 50; i++) {
      checkTumorOnly("A:C", "A:C", contamination, new int[] {100 - i, i, 0, 0});
    }
  }
//  public void testTumorOnly60() {
//    rangeCheck(0.6, 13, 31);
//  }

  public void testTumorOnly70() {
    checkTumorOnly("A:A", "A:C", 0.7, new int[] {85, 15, 0, 0});
  }

  public void testTumorOnly80() {
    //Too close to read error
    checkTumorOnly("A:A", "A:A", 0.7, new int[] {90, 10, 0, 0});
    // Adding more reads reveals the tumor
    checkTumorOnly("A:A", "A:C", 0.7, new int[] {900, 100, 0, 0});
  }
}
