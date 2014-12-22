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
package com.rtg.variant.bayes.multisample.forwardbackward;


import java.util.ArrayList;
import java.util.List;
import java.util.Set;

import com.rtg.relation.Family;
import com.rtg.relation.GenomeRelationships;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.test.NanoRegression;
import com.rtg.variant.GenomePriorParams;
import com.rtg.variant.GenomePriorParamsBuilder;
import com.rtg.variant.VariantOutputLevel;
import com.rtg.variant.VariantParams;
import com.rtg.variant.VariantSample;
import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.ModelInterface;
import com.rtg.variant.bayes.multisample.HaploidDiploidHypotheses;
import com.rtg.variant.bayes.multisample.HypothesisScore;
import com.rtg.variant.bayes.multisample.HypothesisScores;
import com.rtg.variant.bayes.multisample.MultisampleJointCaller;
import com.rtg.variant.bayes.multisample.PriorContainer;
import com.rtg.variant.bayes.multisample.family.AbstractFamilyPosteriorTest;
import com.rtg.variant.bayes.multisample.family.FamilyCallerTest;
import com.rtg.variant.bayes.snp.DescriptionCommon;
import com.rtg.variant.bayes.snp.DescriptionNone;
import com.rtg.variant.bayes.snp.HypothesesNone;
import com.rtg.variant.bayes.snp.HypothesesPrior;
import com.rtg.variant.util.arithmetic.LogApproximatePossibility;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;
import com.rtg.variant.util.arithmetic.SimplePossibility;

/**
 */
public class FamilyCallerFBTest extends FamilyCallerTest {

  @Override
  public void setUp() {
    super.setUp();
    mNano = new NanoRegression(FamilyCallerFBTest.class);
  }

  @Override
  public void tearDown() throws Exception {
    super.tearDown();
  }

  @Override
  protected MultisampleJointCaller getFamilyCaller(Family family, VariantParams vParams) {
    return new FamilyCallerFB(vParams, family) {
      @Override
      protected <D extends Description, T extends HypothesesPrior<D>> ComparisonResult makeSamples(List<ModelInterface<?>> models, HaploidDiploidHypotheses<T> hypotheses) {
        final PriorContainer<T> priorCon = new PriorContainer<>(hypotheses, makeInitialBs(models));
        final HypothesisScores scores = getBestScores(models, priorCon);

        final VariantSample[] samples = new VariantSample[models.size()];
        for (int i = 0; i < samples.length; i++) {
          final HypothesisScore score = scores.getScores()[i];
          if (score != null && score.hypothesis() != -1) {
            final T childHypotheses = hypotheses.get(models.get(i));
            samples[i] = createSample(childHypotheses, score, models.get(i), mParams);
          }
        }

        if (!scores.isInteresting() && mParams.callLevel() != VariantOutputLevel.ALL) {
          return null;
        }

        return new ComparisonResult(scores.isInteresting(), samples, scores.getNonIdentityPosterior());
      }
    };
  }


  //note all resources for this test should be the same as in familycallertest except for the qual field, check using:
  //for f in test/com/rtg/variant/bayes/multisample/family/resources/familycaller*; do diff $f test/com/rtg/variant/bayes/multisample/forwardbackward/resources/familycallerfb-${f##*-}; done
  //and your eyes
  @Override
  protected String nanoPrefix() {
    return "familycallerfb-";
  }

  public void testDiploidDiploidDiploidDenovo() throws InvalidParamsException {
    final GenomePriorParams priors = new GenomePriorParamsBuilder().denovoRef(0.01).denovoNonRef(0.00001).create();
//    final GenomePriorParams priors = new GenomePriorParamsBuilder().denovoRef(0.99).denovoNonRef(0.99).create();
    final PossibilityArithmetic arith = SimplePossibility.SINGLETON;
    final Description desc = new DescriptionCommon("A", "C");
    final HypothesesPrior<Description> none = new HypothesesNone<>(DescriptionNone.SINGLETON, arith, 0);
    final HypothesesPrior<Description> haploid = new AbstractFamilyPosteriorTest.MockHyp(desc, arith, true, new double[] {0.7, 0.3});
    final HypothesesPrior<Description> diploid = new AbstractFamilyPosteriorTest.MockHyp(desc, arith, false, new double[] {0.20, 0.70, 0.10});
    final HaploidDiploidHypotheses<HypothesesPrior<Description>> hdh = new HaploidDiploidHypotheses<>(none, haploid, diploid, false, null);
    final List<ModelInterface<?>> models = new ArrayList<>();
    models.add(new AbstractFamilyPosteriorTest.MockModel(diploid, new double[] {0.1, 0.2, 0.7}));
    models.add(new AbstractFamilyPosteriorTest.MockModel(diploid, new double[] {0.60, 0.25, 0.15}));
    models.add(new AbstractFamilyPosteriorTest.MockModel(diploid, new double[] {0.65, 0.15, 0.20}));
    final Family family = AbstractFamilyPosteriorTest.makeFamilySex(GenomeRelationships.SEX_MALE);
    final FamilyCallerFB ffb = new FamilyCallerFB(VariantParams.builder().genomePriors(priors).create(), family);
    final HypothesisScore[] scores = ffb.getBestScores(models, new PriorContainer<>(hdh, ffb.makeInitialBs(models))).getScores();

    AbstractFamilyPosteriorTest.checkBest(scores[Family.FATHER_INDEX], 1, -0.1286, 1.7263);
    AbstractFamilyPosteriorTest.checkBest(scores[Family.MOTHER_INDEX], 0, 0.1474, -0.1474);
    AbstractFamilyPosteriorTest.checkBest(scores[Family.FIRST_CHILD_INDEX], 2, -0.3677, 0.9046);
  }


  public void testDiploidDiploidDiploidDenovoSon() throws InvalidParamsException {
    final GenomePriorParams priors = new GenomePriorParamsBuilder().denovoRef(2.3E-3).denovoNonRef(2.3E-6).create();
//    final GenomePriorParams priors = new GenomePriorParamsBuilder().denovoRef(0).denovoNonRef(0).create();
    final PossibilityArithmetic arith = SimplePossibility.SINGLETON;
    final Description desc = new DescriptionCommon("A", "C");
    final HypothesesPrior<Description> none = new HypothesesNone<>(DescriptionNone.SINGLETON, arith, 0);
    final HypothesesPrior<Description> haploid = new AbstractFamilyPosteriorTest.MockHyp(desc, arith, true, new double[] {0.7, 0.3});
    final HypothesesPrior<Description> diploid = new AbstractFamilyPosteriorTest.MockHyp(desc, arith, false, new double[] {0.80, 0.10, 0.10});
    final HaploidDiploidHypotheses<HypothesesPrior<Description>> hdh = new HaploidDiploidHypotheses<>(none, haploid, diploid, false, null);
    final List<ModelInterface<?>> models = new ArrayList<>();
    models.add(new AbstractFamilyPosteriorTest.MockModel(diploid, new double[] {0.999, 0.0001, 0.0009}));
    models.add(new AbstractFamilyPosteriorTest.MockModel(diploid, new double[] {0.998, 0.00015, 0.00185}));
    models.add(new AbstractFamilyPosteriorTest.MockModel(diploid, new double[] {0.00005, 0.00005, 0.9999}));
    final Family family = AbstractFamilyPosteriorTest.makeFamilySex(GenomeRelationships.SEX_MALE);
    final FamilyCallerFB ffb = new FamilyCallerFB(VariantParams.builder().genomePriors(priors).create(), family);
    final HypothesisScore[] scores = ffb.getBestScores(models, new PriorContainer<>(hdh, ffb.makeInitialBs(models))).getScores();

    AbstractFamilyPosteriorTest.checkBest(scores[Family.FATHER_INDEX], 0, 2.9647, -2.9647);
    AbstractFamilyPosteriorTest.checkBest(scores[Family.MOTHER_INDEX], 0, 2.2431, -2.2431);
    AbstractFamilyPosteriorTest.checkBest(scores[Family.FIRST_CHILD_INDEX], 2, 3.2981, 3.2981);
    assertEquals(1.5122, scores[Family.FIRST_CHILD_INDEX].getDeNovoPosterior(), 1e-4);
  }

  public void testHaploidNoneHaploid() throws InvalidParamsException {
    final GenomePriorParams priors = new GenomePriorParamsBuilder().denovoRef(0.01).denovoNonRef(0.00001).create();
    final PossibilityArithmetic arith = SimplePossibility.SINGLETON;
    final Description desc = new DescriptionCommon("A", "C");
    final HypothesesPrior<Description> none = new HypothesesNone<>(DescriptionNone.SINGLETON, arith, 0);
    final HypothesesPrior<Description> haploid = new AbstractFamilyPosteriorTest.MockHyp(desc, arith, true, new double[] {0.2, 0.8});
    final HaploidDiploidHypotheses<HypothesesPrior<Description>> hdh = new HaploidDiploidHypotheses<>(none, haploid, null, false, null);
    final List<ModelInterface<?>> models = new ArrayList<>();
    models.add(new AbstractFamilyPosteriorTest.MockModel(haploid, new double[] {0.25, 0.75}));
    models.add(new AbstractFamilyPosteriorTest.MockModel(none));
    models.add(new AbstractFamilyPosteriorTest.MockModel(haploid, new double[] {0.3, 0.7}));
    final Family family = AbstractFamilyPosteriorTest.makeFamilySex(GenomeRelationships.SEX_MALE);
    final FamilyCallerFB ffb = new FamilyCallerFB(VariantParams.builder().genomePriors(priors).create(), family);
    final HypothesisScore[] scores = ffb.getBestScores(models, new PriorContainer<>(hdh, ffb.makeInitialBs(models))).getScores();

    // Not verified by external calculation. just being needlessly more specific than a NaN/Infinity check
    AbstractFamilyPosteriorTest.checkBest(scores[Family.FATHER_INDEX], 1, 3.3245, 3.3244);
    assertNull(scores[Family.MOTHER_INDEX]);
    AbstractFamilyPosteriorTest.checkBest(scores[Family.FIRST_CHILD_INDEX], 1, 3.3324, 3.3324);
  }
  /*
      Fathers left, Mothers right


          A---B
            |
        C---D
          |
          E-----F
            | |
            G H

      Attempt to test pedigree calling on Y chromosome
   */
  public void testPedigreeHaploidNoneHaploid() throws InvalidParamsException {
    final GenomePriorParams priors = new GenomePriorParamsBuilder().denovoRef(0.01).denovoNonRef(0.00001).create();

    final GenomeRelationships pedigree = new GenomeRelationships();
    pedigree.addGenome("A", GenomeRelationships.SEX_MALE);
    pedigree.addGenome("B", GenomeRelationships.SEX_FEMALE);
    pedigree.addGenome("C", GenomeRelationships.SEX_MALE);
    pedigree.addGenome("D", GenomeRelationships.SEX_FEMALE);
    pedigree.addGenome("E", GenomeRelationships.SEX_MALE);
    pedigree.addGenome("F", GenomeRelationships.SEX_FEMALE);
    pedigree.addGenome("G", GenomeRelationships.SEX_MALE);
    pedigree.addGenome("H", GenomeRelationships.SEX_FEMALE);

    pedigree.addParentChild("A", "D");
    pedigree.addParentChild("B", "D");
    pedigree.addParentChild("C", "E");
    pedigree.addParentChild("D", "E");
    pedigree.addParentChild("E", "G");
    pedigree.addParentChild("F", "G");
    pedigree.addParentChild("E", "H");
    pedigree.addParentChild("F", "H");

    final PossibilityArithmetic arith = LogApproximatePossibility.SINGLETON;
    final Description desc = new DescriptionCommon("A", "C");

    final HypothesesPrior<Description> none = new HypothesesNone<>(DescriptionNone.SINGLETON, arith, 0);
    final HypothesesPrior<Description> haploid = new AbstractFamilyPosteriorTest.MockHyp(desc, arith, true, new double[] {0.2, 0.8});
    final HaploidDiploidHypotheses<HypothesesPrior<Description>> hdh = new HaploidDiploidHypotheses<>(none, haploid, null, false, null);

    final List<ModelInterface<?>> models = new ArrayList<>();
    models.add(new AbstractFamilyPosteriorTest.MockModel(haploid, new double[]{0.25, 0.75}));
    models.add(new AbstractFamilyPosteriorTest.MockModel(none));
    models.add(new AbstractFamilyPosteriorTest.MockModel(haploid, new double[]{0.3, 0.7}));
    models.add(new AbstractFamilyPosteriorTest.MockModel(none));
    models.add(new AbstractFamilyPosteriorTest.MockModel(haploid, new double[]{0.3, 0.7}));
    models.add(new AbstractFamilyPosteriorTest.MockModel(none));
    models.add(new AbstractFamilyPosteriorTest.MockModel(haploid, new double[]{0.25, 0.75}));
    models.add(new AbstractFamilyPosteriorTest.MockModel(none));

    final Set<Family> families = Family.getFamilies(pedigree, false, null);
    final List<String> calledGenomes = new ArrayList<>();
    // Ensure all family members are in called genomes
    for (Family family : families) {
      for (String member : family.getMembers()) {
        if (!calledGenomes.contains(member)) {
          calledGenomes.add(member);
        }
      }
    }
    for (Family family : families) {
      family.setSampleIds(calledGenomes);
    }
    final FamilyCallerFB ffb = new FamilyCallerFB(VariantParams.builder().genomePriors(priors).create()
        , families.toArray(new Family[families.size()])
        );
    final HypothesisScores bestScores = ffb.getBestScores(models, new PriorContainer<>(hdh, ffb.makeInitialBs(models)));
    assertFalse(Double.isNaN(bestScores.getNonIdentityPosterior()));
    assertFalse(Double.isInfinite(bestScores.getNonIdentityPosterior()));
    final HypothesisScore[] scores = bestScores.getScores();
    checkScores(scores);

  }

  // Simple NaN / Infinity check
  public static void checkScores(HypothesisScore[] scores) {
    for (HypothesisScore score : scores) {
      // When HypothesesNone Will be null if coming from family caller, -1 when coming from pop caller
      if (score != null && score.hypothesis() != -1) {
        assertFalse(Double.isNaN(score.posterior()));
        assertFalse(Double.isInfinite(score.posterior()));
        assertFalse(Double.isNaN(score.nonIdentityPosterior()));
        assertFalse(Double.isInfinite(score.nonIdentityPosterior()));
      }
    }
  }

}
