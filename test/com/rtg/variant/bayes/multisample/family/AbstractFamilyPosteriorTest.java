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

package com.rtg.variant.bayes.multisample.family;


import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import com.rtg.relation.Family;
import com.rtg.relation.GenomeRelationships;
import com.rtg.util.InvalidParamsException;
import com.rtg.variant.GenomePriorParams;
import com.rtg.variant.GenomePriorParamsBuilder;
import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.Hypotheses;
import com.rtg.variant.bayes.Model;
import com.rtg.variant.bayes.ModelInterface;
import com.rtg.variant.bayes.NoAlleleBalance;
import com.rtg.variant.bayes.multisample.HaploidDiploidHypotheses;
import com.rtg.variant.bayes.multisample.HypothesisScore;
import com.rtg.variant.bayes.multisample.forwardbackward.FamilyPosteriorFBTest;
import com.rtg.variant.bayes.snp.DescriptionCommon;
import com.rtg.variant.bayes.snp.DescriptionNone;
import com.rtg.variant.bayes.snp.HypothesesMock;
import com.rtg.variant.bayes.snp.HypothesesNone;
import com.rtg.variant.bayes.snp.HypothesesPrior;
import com.rtg.variant.bayes.snp.StatisticsSnp;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;
import com.rtg.variant.util.arithmetic.SimplePossibility;

import junit.framework.TestCase;

/**
 */
public abstract class AbstractFamilyPosteriorTest extends TestCase {
  protected static final String MOTHER = "m";
  protected static final String FATHER = "f";

  private HashMap<String, Double> mExpectedNonIdentity;
  @Override
  public void setUp() throws Exception {
    super.setUp();
    mExpectedNonIdentity = new HashMap<>();
    mExpectedNonIdentity.put("testHaploidNoneHaploid", -0.5596);
    mExpectedNonIdentity.put("testNoneHaploidHaploid", -0.5596);
    mExpectedNonIdentity.put("testHaploidDiploidHaploid", 0.4564);
    mExpectedNonIdentity.put("testHaploidDiploidHaploidAvian", 0.4564);
    mExpectedNonIdentity.put("testHaploidDiploidDiploid", -0.0344);
    mExpectedNonIdentity.put("testHaploidDiploidDiploidAvian", -0.0344);
    mExpectedNonIdentity.put("testHaploidHaploidHaploid", -0.1028); // <- not externally validated
    mExpectedNonIdentity.put("testDiploidDiploidDiploid", 2.2127);
    mExpectedNonIdentity.put("testDiploidDiploidDiploidSon", 2.2127);
    mExpectedNonIdentity.put("testDiploidDiploidDiploidDiploid", 1.8427);
    mExpectedNonIdentity.put("testDiploidDiploidDiploidDiploidDiploid", 2.1063);
    mExpectedNonIdentity.put("testHaploidNoneHaploidNoneHaploid", -1.9459);
    mExpectedNonIdentity.put("testHaploidDiploidHaploidDiploidHaploid", 0.6336);
    mExpectedNonIdentity.put("testDiploidDiploidDiploidDenovo", 2.2129);
    mExpectedNonIdentity.put("testDiploidDiploidDiploidDenovoSon", 3.2983);
    mExpectedNonIdentity.put("testDiploidDiploidDiploidDenovoSonNoPrior", 1.4035);
  }

  @Override
  protected void tearDown() throws Exception {
    super.tearDown();
    mExpectedNonIdentity = null;
  }

  protected double getExpectedNonIdentity(String method) {
    return mExpectedNonIdentity.get(method);
  }


  // Following tests are tied back to spreadsheet familyposteriortest.xlsx and exercise the various sex chromosomes and
  // assume priors are supplied explicitly

  /**
   * Mock hypotheses class
   */
  public static class MockHyp extends HypothesesMock<Description> {
    public MockHyp(Description desc, PossibilityArithmetic arith, boolean haploid, double[] priors) {
      super(desc, arith, haploid, 0, priors);
    }
    @Override
    protected void initPriors(double[] priors) {
      //    do nothing
    }
  }

  /**
   * Mock model class
   */
  public static class MockModel extends Model<Description> {
    public MockModel(final Hypotheses<Description> hyp, final double[] posteriors) {
      super(hyp, new StatisticsSnp(hyp.description()), new NoAlleleBalance());
      assert posteriors.length == mPosteriors.length;
      for (int i = 0; i < posteriors.length; ++i) {
        mPosteriors[i] = arithmetic().prob2Poss(posteriors[i]);
      }
    }

    public MockModel(final Hypotheses<Description> hyp) {
      super(hyp, new StatisticsSnp(hyp.description()), new NoAlleleBalance());
      for (int i = 0; i < mPosteriors.length; ++i) {
        mPosteriors[i] = hyp.arithmetic().one();
      }
    }


  }

  public static Family makeFamilySex(String...children) {
    final GenomeRelationships pedigree = new GenomeRelationships();
    pedigree.addGenome(FATHER, GenomeRelationships.SEX_MALE);
    pedigree.addGenome(MOTHER, GenomeRelationships.SEX_FEMALE);
    final String[] childNames = new String[children.length];
    for (int i = 0; i < children.length; ++i) {
      final String child = "c" + i;
      childNames[i] = child;
      pedigree.addGenome(child, children[i]);
      pedigree.addParentChild(FATHER, child);
      pedigree.addParentChild(MOTHER, child);
    }
    //System.err.println(pedigree);
    return new Family(pedigree, FATHER, MOTHER, childNames);
  }


  protected AbstractFamilyPosterior getFamilyPosterior(final GenomePriorParams priors, final HaploidDiploidHypotheses<HypothesesPrior<Description>> hdh, final List<ModelInterface<?>> models, final Family family) {
    for (ModelInterface<?> model : models) {
      model.freeze();
    }
    return new FamilyPosterior(family, priors, models, hdh);
  }

  public static void checkBest(final HypothesisScore best, final int call, final double score, final double auxScore) {
    //System.err.println(best);
    assertEquals(call, best.hypothesis());
    assertEquals(score, best.posterior(), 1.5e-4);
    assertEquals(auxScore, best.nonIdentityPosterior(), 1.5e-4);
  }

  //Template for start of tests
  //  final GenomePriorParams priors = new GenomePriorParamsBuilder().create();
  //  final PossibilityArithmetic arith = SimplePossibility.SINGLETON;
  //  final Description desc = new DescriptionCommon("A", "C");
  //  final Hypotheses<Description> none = new HypothesesNone(arith, priors, 0);
  //  final Hypotheses<Description> haploid = new MockHyp(desc, arith, true, new double[] {0.0, 0.0}); ///TODO priors
  //  final Hypotheses<Description> diploid = new MockHyp(desc, arith, false, new double[] {0.0, 0.0, 0.0}); ///TODO priors
  //  final HaploidDiploidHypotheses<Hypotheses<Description>> hdh = new HaploidDiploidHypotheses<>(haploid, diploid, false);
  //  final List<ModelInterface<?>> models = new ArrayList<>();
  //  models.add(new MockModel(diploid, new double[] {0.0, 0.0, 0.0}));  ///TODO priors
  //  final Family family = makeFamilySex(GenomeRelationships.SEX_MALE);
  //  final FamilyPosterior fp = new FamilyPosterior(family, priors, models, hdh);

  //human Y chromosome
  public void testHaploidNoneHaploid() throws InvalidParamsException {
    final GenomePriorParams priors = getGenomePriorParams();
    final PossibilityArithmetic arith = SimplePossibility.SINGLETON;
    final Description desc = new DescriptionCommon("A", "C");
    final HypothesesPrior<Description> none = new HypothesesNone<>(DescriptionNone.SINGLETON, arith, 0);
    final HypothesesPrior<Description> haploid = new MockHyp(desc, arith, true, new double[] {0.2, 0.8});
    final HaploidDiploidHypotheses<HypothesesPrior<Description>> hdh = new HaploidDiploidHypotheses<>(none, haploid, null, false, null);
    final List<ModelInterface<?>> models = new ArrayList<>();
    models.add(new MockModel(haploid, new double[] {0.75, 0.25}));
    models.add(new MockModel(none));
    models.add(new MockModel(haploid, new double[] {0.7, 0.3}));
    final Family family = makeFamilySex(GenomeRelationships.SEX_MALE);
    final AbstractFamilyPosterior fp = getFamilyPosterior(priors, hdh, models, family);
    checkBest(fp.bestFather(), 0, 0.5596, -0.5596);
    assertNull(fp.bestMother());
    checkBest(fp.bestChild(0), 0, 0.5596, -0.5596);
    assertEquals(1, fp.numberChildren());
    assertFalse(fp.isInteresting());
    assertEquals(getExpectedNonIdentity("testHaploidNoneHaploid"), fp.getNonIdentityPosterior(), 1e-4);
  }
  //human mitochondrial
  public void testHaploidHaploidHaploid() throws InvalidParamsException {
    final AbstractFamilyPosterior fp = getHaploidHaploidHaploid(false);
    checkBest(fp.bestFather(), 0, 1.5581, -1.5581); // TODO compute this externally. This is currently just testing that it can make calls
    checkBest(fp.bestMother(), 0, 0.5596, -0.5596);
    checkBest(fp.bestChild(0), 0, 0.5596, -0.5596);
    assertEquals(1, fp.numberChildren());
    assertFalse(fp.isInteresting());
    assertEquals(getExpectedNonIdentity("testHaploidHaploidHaploid"), fp.getNonIdentityPosterior(), 1e-4);
  }
  //human mitochondrial
  public void testHaploidHaploidHaploidNonDefault() throws InvalidParamsException {
    final boolean isDefault = true;
    // Check that this is runs ok with default priors.
    getHaploidHaploidHaploid(isDefault);
  }

  private AbstractFamilyPosterior getHaploidHaploidHaploid(boolean aDefault) {
    return getHaploidHaploidHaploid(aDefault, new double[]{0.95, 0.05}, new double[]{0.75, 0.25});
  }
    private AbstractFamilyPosterior getHaploidHaploidHaploid(boolean aDefault, double[] father, double[] mother) {
    final GenomePriorParams priors = getGenomePriorParams();
    final PossibilityArithmetic arith = SimplePossibility.SINGLETON;
    final Description desc = new DescriptionCommon("A", "C");
    final HypothesesPrior<Description> haploid = new MockHyp(desc, arith, true, new double[] {0.2, 0.8});
    final HypothesesPrior<Description> none = new HypothesesNone<>(DescriptionNone.SINGLETON, arith, 0);
    final HaploidDiploidHypotheses<HypothesesPrior<Description>> hdh = new HaploidDiploidHypotheses<>(none, haploid, null, aDefault, null);
    final List<ModelInterface<?>> models = new ArrayList<>();
    models.add(new MockModel(haploid, father));
    models.add(new MockModel(haploid, mother));
    models.add(new MockModel(haploid, new double[] {0.7, 0.3}));
    final Family family = makeFamilySex(GenomeRelationships.SEX_MALE);
    return getFamilyPosterior(priors, hdh, models, family);
  }
  //human mitochondrial
  public void testHaploidHaploidHaploidFatherIrrelevant() throws InvalidParamsException {
    final AbstractFamilyPosterior fp = getHaploidHaploidHaploid(false, new double[]{0.1, 0.9}, new double[]{0.75, 0.25});
    final AbstractFamilyPosterior fp2 = getHaploidHaploidHaploid(false, new double[]{0.9, 0.1}, new double[]{0.75, 0.25});
    assertEquals(fp.bestChild(0).posterior(), fp2.bestChild(0).posterior(), 0.001);
  }
  //human mitochondrial
  public void testHaploidHaploidHaploidMotherRelevant() throws InvalidParamsException {
    final AbstractFamilyPosterior fp = getHaploidHaploidHaploid(false, new double[]{0.1, 0.9}, new double[]{0.75, 0.25});
    final AbstractFamilyPosterior fp2 = getHaploidHaploidHaploid(false, new double[]{0.9, 0.1}, new double[]{0.99, 0.01});
    final double posteriorLow = fp.bestChild(0).posterior();
    final double posteriorHigh = fp2.bestChild(0).posterior();
    final String err = posteriorLow + " < " + posteriorHigh;
    assertTrue(err, posteriorLow < posteriorHigh - 3);
  }

  //human Y chromosome
  public void testHaploidNoneHaploidNoneHaploid() throws InvalidParamsException {
    final GenomePriorParams priors = getGenomePriorParams();
    final PossibilityArithmetic arith = SimplePossibility.SINGLETON;
    final Description desc = new DescriptionCommon("A", "C");
    final HypothesesPrior<Description> none = new HypothesesNone<>(DescriptionNone.SINGLETON, arith, 0);
    final HypothesesPrior<Description> haploid = new MockHyp(desc, arith, true, new double[] {0.2, 0.8});
    final HaploidDiploidHypotheses<HypothesesPrior<Description>> hdh = new HaploidDiploidHypotheses<>(none, haploid, null, false, null);
    final List<ModelInterface<?>> models = new ArrayList<>();
    models.add(new MockModel(haploid, new double[] {0.75, 0.25}));
    models.add(new MockModel(none));
    models.add(new MockModel(haploid, new double[] {0.7, 0.3}));
    models.add(new MockModel(none));
    models.add(new MockModel(haploid, new double[] {0.8, 0.2}));
    final Family family = makeFamilySex(GenomeRelationships.SEX_MALE, GenomeRelationships.SEX_FEMALE, GenomeRelationships.SEX_MALE);
    final AbstractFamilyPosterior fp = getFamilyPosterior(priors, hdh, models, family);
    checkBest(fp.bestFather(), 0, 1.9459, -1.9459);
    assertNull(fp.bestMother());
    checkBest(fp.bestChild(0), 0, 1.9459, -1.9459);
    assertNull(fp.bestChild(1));
    checkBest(fp.bestChild(2), 0, 1.9459, -1.9459);
    assertEquals(3, fp.numberChildren());
    assertFalse(fp.isInteresting());
    assertEquals(getExpectedNonIdentity("testHaploidNoneHaploidNoneHaploid"), fp.getNonIdentityPosterior(), 1e-4);
  }





  //avian W chromosome
  public void testNoneHaploidHaploid() throws InvalidParamsException {
    final GenomePriorParams priors = getGenomePriorParams();
    final PossibilityArithmetic arith = SimplePossibility.SINGLETON;
    final Description desc = new DescriptionCommon("A", "C");
    final HypothesesPrior<Description> none = new HypothesesNone<>(DescriptionNone.SINGLETON, arith, 0);
    final HypothesesPrior<Description> haploid = new MockHyp(desc, arith, true, new double[] {0.2, 0.8});
    final HaploidDiploidHypotheses<HypothesesPrior<Description>> hdh = new HaploidDiploidHypotheses<>(none, haploid, null, false, null);
    final List<ModelInterface<?>> models = new ArrayList<>();
    models.add(new MockModel(none));
    models.add(new MockModel(haploid, new double[] {0.75, 0.25}));
    models.add(new MockModel(haploid, new double[] {0.7, 0.3}));
    final Family family = makeFamilySex(GenomeRelationships.SEX_FEMALE);
    final AbstractFamilyPosterior fp = getFamilyPosterior(priors, hdh, models, family);
    checkBest(fp.bestMother(), 0, 0.5596, -0.5596);
    assertNull(fp.bestFather());
    checkBest(fp.bestChild(0), 0, 0.5596, -0.5596);
    assertEquals(1, fp.numberChildren());
    assertFalse(fp.isInteresting());
    assertEquals(getExpectedNonIdentity("testNoneHaploidHaploid"), fp.getNonIdentityPosterior(), 1e-4);
  }

  //human X chromosome with son
  public void testHaploidDiploidHaploid() throws InvalidParamsException {
    final GenomePriorParams priors = getGenomePriorParams();
    final PossibilityArithmetic arith = SimplePossibility.SINGLETON;
    final Description desc = new DescriptionCommon("A", "C");
    final HypothesesPrior<Description> none = new HypothesesNone<>(DescriptionNone.SINGLETON, arith, 0);
    final HypothesesPrior<Description> haploid = new MockHyp(desc, arith, true, new double[] {0.7, 0.3});
    final HypothesesPrior<Description> diploid = new MockHyp(desc, arith, false, new double[] {0.20, 0.70, 0.10});
    final HaploidDiploidHypotheses<HypothesesPrior<Description>> hdh = new HaploidDiploidHypotheses<>(none, haploid, diploid, false, null);
    final List<ModelInterface<?>> models = new ArrayList<>();
    models.add(new MockModel(haploid, new double[] {0.4, 0.6}));
    models.add(new MockModel(diploid, new double[] {0.60, 0.25, 0.15}));
    models.add(new MockModel(haploid, new double[] {0.75, 0.25}));
    final Family family = makeFamilySex(GenomeRelationships.SEX_MALE);
    final AbstractFamilyPosterior fp = getFamilyPosterior(priors, hdh, models, family);
    checkBest(fp.bestFather(), 0, 0.4418, -0.4418);
    checkBest(fp.bestMother(), 0, 0.5631, -0.5631/*??*/);
    checkBest(fp.bestChild(0), 0, 0.7400, -0.7400);
    assertEquals(1, fp.numberChildren());
    assertFalse(fp.isInteresting());
    assertEquals(getExpectedNonIdentity("testHaploidDiploidHaploid"), fp.getNonIdentityPosterior(), 1e-4);
  }

    //human X chromosome with son and daughter and son again
  public void testHaploidDiploidHaploidDiploidHaploid() throws InvalidParamsException {
    final GenomePriorParams priors = getGenomePriorParams();
    final PossibilityArithmetic arith = SimplePossibility.SINGLETON;
    final Description desc = new DescriptionCommon("A", "C");
    final HypothesesPrior<Description> none = new HypothesesNone<>(DescriptionNone.SINGLETON, arith, 0);
    final HypothesesPrior<Description> haploid = new MockHyp(desc, arith, true, new double[] {0.7, 0.3});
    final HypothesesPrior<Description> diploid = new MockHyp(desc, arith, false, new double[] {0.20, 0.70, 0.10});
    final HaploidDiploidHypotheses<HypothesesPrior<Description>> hdh = new HaploidDiploidHypotheses<>(none, haploid, diploid, false, null);
    final List<ModelInterface<?>> models = new ArrayList<>();
    models.add(new MockModel(haploid, new double[] {0.4, 0.6}));
    models.add(new MockModel(diploid, new double[] {0.60, 0.25, 0.15}));
    models.add(new MockModel(haploid, new double[] {0.75, 0.25}));
    models.add(new MockModel(diploid, new double[] {0.45, 0.20, 0.35}));
    models.add(new MockModel(haploid, new double[] {0.30, 0.70}));
    final Family family = makeFamilySex(GenomeRelationships.SEX_MALE, GenomeRelationships.SEX_FEMALE, GenomeRelationships.SEX_MALE);
    final AbstractFamilyPosterior fp = getFamilyPosterior(priors, hdh, models, family);
    checkBest(fp.bestFather(), 0, 0.8254, -0.8254);
    checkBest(fp.bestMother(), 0, 0.0803, -0.0803);
    checkBest(fp.bestChild(0), 0, 0.2670, -0.2670);
    checkBest(fp.bestChild(1), 2, 0.0399, 0.5289);
    checkBest(fp.bestChild(2), 0, 0.1546, -0.1546);
    assertEquals(3, fp.numberChildren());
    assertTrue(fp.isInteresting());
    assertEquals(getExpectedNonIdentity("testHaploidDiploidHaploidDiploidHaploid"), fp.getNonIdentityPosterior(), 1e-4);
  }


  //avian Z chromosome with daughter
  public void testHaploidDiploidHaploidAvian() throws InvalidParamsException {
    final GenomePriorParams priors = getGenomePriorParams();
    final PossibilityArithmetic arith = SimplePossibility.SINGLETON;
    final Description desc = new DescriptionCommon("A", "C");
    final HypothesesPrior<Description> none = new HypothesesNone<>(DescriptionNone.SINGLETON, arith, 0);
    final HypothesesPrior<Description> haploid = new MockHyp(desc, arith, true, new double[] {0.7, 0.3});
    final HypothesesPrior<Description> diploid = new MockHyp(desc, arith, false, new double[] {0.20, 0.70, 0.10});
    final HaploidDiploidHypotheses<HypothesesPrior<Description>> hdh = new HaploidDiploidHypotheses<>(none, haploid, diploid, false, null);
    final List<ModelInterface<?>> models = new ArrayList<>();
    models.add(new MockModel(diploid, new double[] {0.60, 0.25, 0.15}));
    models.add(new MockModel(haploid, new double[] {0.4, 0.6}));
    models.add(new MockModel(haploid, new double[] {0.75, 0.25}));
    final Family family = makeFamilySex(GenomeRelationships.SEX_FEMALE);
    final AbstractFamilyPosterior fp = getFamilyPosterior(priors, hdh, models, family);
    checkBest(fp.bestFather(), 0, 0.5631, -0.5631/*??*/);
    checkBest(fp.bestMother(), 0, 0.4418, -0.4418);
    checkBest(fp.bestChild(0), 0, 0.7400, -0.7400);
    assertEquals(1, fp.numberChildren());
    assertFalse(fp.isInteresting());
    assertEquals(getExpectedNonIdentity("testHaploidDiploidHaploidAvian"), fp.getNonIdentityPosterior(), 1e-4);
  }

  //human X chromosome with daughter
  public void testHaploidDiploidDiploid() throws InvalidParamsException {
    final GenomePriorParams priors = getGenomePriorParams();
    final PossibilityArithmetic arith = SimplePossibility.SINGLETON;
    final Description desc = new DescriptionCommon("A", "C");
    final HypothesesPrior<Description> none = new HypothesesNone<>(DescriptionNone.SINGLETON, arith, 0);
    final HypothesesPrior<Description> haploid = new MockHyp(desc, arith, true, new double[] {0.7, 0.3});
    final HypothesesPrior<Description> diploid = new MockHyp(desc, arith, false, new double[] {0.20, 0.70, 0.10});
    final HaploidDiploidHypotheses<HypothesesPrior<Description>> hdh = new HaploidDiploidHypotheses<>(none, haploid, diploid, false, null);
    final List<ModelInterface<?>> models = new ArrayList<>();
    models.add(new MockModel(haploid, new double[] {0.4, 0.6}));
    models.add(new MockModel(diploid, new double[] {0.60, 0.25, 0.15}));
    models.add(new MockModel(diploid, new double[] {0.65, 0.15, 0.20}));
    final Family family = makeFamilySex(GenomeRelationships.SEX_FEMALE);
    final AbstractFamilyPosterior fp = getFamilyPosterior(priors, hdh, models, family);
    checkBest(fp.bestFather(), 0, 1.2562, -1.2562);
    checkBest(fp.bestMother(), 0, 0.4439, -0.4439);
    checkBest(fp.bestChild(0), 0, 0.1618, -0.1618);
    assertEquals(1, fp.numberChildren());
    assertFalse(fp.isInteresting());
    assertEquals(getExpectedNonIdentity("testHaploidDiploidDiploid"), fp.getNonIdentityPosterior(), 1e-4);
  }

  //avian Z chromosome with son
  public void testHaploidDiploidDiploidAvian() throws InvalidParamsException {
    final GenomePriorParams priors = getGenomePriorParams();
    final PossibilityArithmetic arith = SimplePossibility.SINGLETON;
    final Description desc = new DescriptionCommon("A", "C");
    final HypothesesPrior<Description> none = new HypothesesNone<>(DescriptionNone.SINGLETON, arith, 0);
    final HypothesesPrior<Description> haploid = new MockHyp(desc, arith, true, new double[] {0.7, 0.3});
    final HypothesesPrior<Description> diploid = new MockHyp(desc, arith, false, new double[] {0.20, 0.70, 0.10});
    final HaploidDiploidHypotheses<HypothesesPrior<Description>> hdh = new HaploidDiploidHypotheses<>(none, haploid, diploid, false, null);
    final List<ModelInterface<?>> models = new ArrayList<>();
    models.add(new MockModel(diploid, new double[]{0.60, 0.25, 0.15}));
    models.add(new MockModel(haploid, new double[]{0.4, 0.6}));
    models.add(new MockModel(diploid, new double[]{0.65, 0.15, 0.20}));
    final Family family = makeFamilySex(GenomeRelationships.SEX_MALE);
    final AbstractFamilyPosterior fp = getFamilyPosterior(priors, hdh, models, family);
    checkBest(fp.bestFather(), 0, 0.4439, -0.4439);
    checkBest(fp.bestMother(), 0, 1.2562, -1.2562);
    checkBest(fp.bestChild(0), 0, 0.1618, -0.1618);
    assertEquals(1, fp.numberChildren());
    assertFalse(fp.isInteresting());
    assertEquals(getExpectedNonIdentity("testHaploidDiploidDiploidAvian"), fp.getNonIdentityPosterior(), 1e-4);
  }


    //all autosomes daughter
  public void testDiploidDiploidDiploidDiploid() throws InvalidParamsException {
    final GenomePriorParams priors = getGenomePriorParams();
    final PossibilityArithmetic arith = SimplePossibility.SINGLETON;
    final Description desc = new DescriptionCommon("A", "C");
    final HypothesesPrior<Description> none = new HypothesesNone<>(DescriptionNone.SINGLETON, arith, 0);
    final HypothesesPrior<Description> haploid = new MockHyp(desc, arith, true, new double[] {0.7, 0.3});
    final HypothesesPrior<Description> diploid = new MockHyp(desc, arith, false, new double[] {0.20, 0.70, 0.10});
    final HaploidDiploidHypotheses<HypothesesPrior<Description>> hdh = new HaploidDiploidHypotheses<>(none, haploid, diploid, false, null);
    final List<ModelInterface<?>> models = new ArrayList<>();
    models.add(new MockModel(diploid, new double[] {0.1, 0.2, 0.7}));
    models.add(new MockModel(diploid, new double[] {0.60, 0.25, 0.15}));
    models.add(new MockModel(diploid, new double[] {0.65, 0.15, 0.20}));
    models.add(new MockModel(diploid, new double[] {0.45, 0.20, 0.35}));
    final Family family = makeFamilySex(GenomeRelationships.SEX_FEMALE, GenomeRelationships.SEX_FEMALE);
    final AbstractFamilyPosterior fp = getFamilyPosterior(priors, hdh, models, family);
    checkBest(fp.bestFather(), 2, -0.3487, 1.4217);
    checkBest(fp.bestMother(), 0, 0.5925, -0.5926); // NOTE off by 0.0001
    checkBest(fp.bestChild(0), 2, -0.2886, 0.5397);
    checkBest(fp.bestChild(1), 2, -0.0033, 0.8203 /*should be 0.8204 but FastFamily is slightly off*/);
    assertEquals(2, fp.numberChildren());
    assertTrue(fp.isInteresting());
    assertEquals(getExpectedNonIdentity("testDiploidDiploidDiploidDiploid"), fp.getNonIdentityPosterior(), 1e-4);
  }
      //all autosomes daughter
  public void testDiploidDiploidDiploidDiploidDiploid() throws InvalidParamsException {
    final GenomePriorParams priors = getGenomePriorParams();
    final PossibilityArithmetic arith = SimplePossibility.SINGLETON;
    final Description desc = new DescriptionCommon("A", "C");
    final HypothesesPrior<Description> none = new HypothesesNone<>(DescriptionNone.SINGLETON, arith, 0);
    final HypothesesPrior<Description> haploid = new MockHyp(desc, arith, true, new double[] {0.7, 0.3});
    final HypothesesPrior<Description> diploid = new MockHyp(desc, arith, false, new double[] {0.20, 0.70, 0.10});
    final HaploidDiploidHypotheses<HypothesesPrior<Description>> hdh = new HaploidDiploidHypotheses<>(none, haploid, diploid, false, null);
    final List<ModelInterface<?>> models = new ArrayList<>();
    models.add(new MockModel(diploid, new double[] {0.1, 0.2, 0.7}));
    models.add(new MockModel(diploid, new double[] {0.60, 0.25, 0.15}));
    models.add(new MockModel(diploid, new double[] {0.65, 0.15, 0.20}));
    models.add(new MockModel(diploid, new double[] {0.45, 0.20, 0.35}));
    models.add(new MockModel(diploid, new double[] {0.3, 0.20, 0.5}));
    final Family family = makeFamilySex(GenomeRelationships.SEX_FEMALE, GenomeRelationships.SEX_FEMALE, GenomeRelationships.SEX_FEMALE);
    final AbstractFamilyPosterior fp = getFamilyPosterior(priors, hdh, models, family);
    checkBest(fp.bestFather(), 2, -0.3184, 1.5022);
    checkBest(fp.bestMother(), 0, 0.8742, -0.8742);
    checkBest(fp.bestChild(0), 2, 0.0628, 0.6074);
    checkBest(fp.bestChild(1), 2, 0.3607, 0.9112);
    checkBest(fp.bestChild(2), 2, 0.6592, 1.2261);
    assertEquals(3, fp.numberChildren());
    assertTrue(fp.isInteresting());
    assertEquals(getExpectedNonIdentity("testDiploidDiploidDiploidDiploidDiploid"), fp.getNonIdentityPosterior(), 1e-4);
  }



  //all autosomes daughter
  public void testDiploidDiploidDiploid() throws InvalidParamsException {
    final GenomePriorParams priors = getGenomePriorParams();
    final PossibilityArithmetic arith = SimplePossibility.SINGLETON;
    final Description desc = new DescriptionCommon("A", "C");
    final HypothesesPrior<Description> none = new HypothesesNone<>(DescriptionNone.SINGLETON, arith, 0);
    final HypothesesPrior<Description> haploid = new MockHyp(desc, arith, true, new double[] {0.7, 0.3});
    final HypothesesPrior<Description> diploid = new MockHyp(desc, arith, false, new double[] {0.20, 0.70, 0.10});
    final HaploidDiploidHypotheses<HypothesesPrior<Description>> hdh = new HaploidDiploidHypotheses<>(none, haploid, diploid, false, null);
    final List<ModelInterface<?>> models = new ArrayList<>();
    models.add(new MockModel(diploid, new double[] {0.1, 0.2, 0.7}));
    models.add(new MockModel(diploid, new double[] {0.60, 0.25, 0.15}));
    models.add(new MockModel(diploid, new double[] {0.65, 0.15, 0.20}));
    final Family family = makeFamilySex(GenomeRelationships.SEX_FEMALE);
    final AbstractFamilyPosterior fp = getFamilyPosterior(priors, hdh, models, family);
    checkBest(fp.bestFather(), 1, -0.1283, 1.7273);
    checkBest(fp.bestMother(), 0, 0.1471, -0.1471);
    checkBest(fp.bestChild(0), 2, -0.3681, 0.9044);
    assertEquals(1, fp.numberChildren());
    assertTrue(fp.isInteresting());
    assertEquals(getExpectedNonIdentity("testDiploidDiploidDiploid"), fp.getNonIdentityPosterior(), 1e-4);
  }


  //all autosomes son
  public void testDiploidDiploidDiploidSon() throws InvalidParamsException {
    final GenomePriorParams priors = getGenomePriorParams();
    final PossibilityArithmetic arith = SimplePossibility.SINGLETON;
    final Description desc = new DescriptionCommon("A", "C");
    final HypothesesPrior<Description> none = new HypothesesNone<>(DescriptionNone.SINGLETON, arith, 0);
    final HypothesesPrior<Description> haploid = new MockHyp(desc, arith, true, new double[] {0.7, 0.3});
    final HypothesesPrior<Description> diploid = new MockHyp(desc, arith, false, new double[] {0.20, 0.70, 0.10});
    final HaploidDiploidHypotheses<HypothesesPrior<Description>> hdh = new HaploidDiploidHypotheses<>(none, haploid, diploid, false, null);
    final List<ModelInterface<?>> models = new ArrayList<>();
    models.add(new MockModel(diploid, new double[] {0.1, 0.2, 0.7}));
    models.add(new MockModel(diploid, new double[] {0.60, 0.25, 0.15}));
    models.add(new MockModel(diploid, new double[] {0.65, 0.15, 0.20}));
    final Family family = makeFamilySex(GenomeRelationships.SEX_MALE);
    final AbstractFamilyPosterior fp = getFamilyPosterior(priors, hdh, models, family);
    checkBest(fp.bestFather(), 1, -0.1283, 1.7273);
    checkBest(fp.bestMother(), 0, 0.1471, -0.1471);
    checkBest(fp.bestChild(0), 2, -0.3681, 0.9044);
    assertEquals(1, fp.numberChildren());
    assertTrue(fp.isInteresting());
    assertEquals(getExpectedNonIdentity("testDiploidDiploidDiploidSon"), fp.getNonIdentityPosterior(), 1e-4);
  }

  //all autosomes with denovo priors
  public void testDiploidDiploidDiploidDenovo() throws InvalidParamsException {
    final GenomePriorParams priors = new GenomePriorParamsBuilder().denovoRef(0.01).denovoNonRef(0.00001).create();
//    final GenomePriorParams priors = new GenomePriorParamsBuilder().denovoRef(0.99).denovoNonRef(0.99).create();
    final PossibilityArithmetic arith = SimplePossibility.SINGLETON;
    final Description desc = new DescriptionCommon("A", "C");
    final HypothesesPrior<Description> none = new HypothesesNone<>(DescriptionNone.SINGLETON, arith, 0);
    final HypothesesPrior<Description> haploid = new MockHyp(desc, arith, true, new double[] {0.7, 0.3});
    final HypothesesPrior<Description> diploid = new MockHyp(desc, arith, false, new double[] {0.20, 0.70, 0.10});
    final HaploidDiploidHypotheses<HypothesesPrior<Description>> hdh = new HaploidDiploidHypotheses<>(none, haploid, diploid, false, null);
    final List<ModelInterface<?>> models = new ArrayList<>();
    models.add(new MockModel(diploid, new double[] {0.1, 0.2, 0.7}));
    models.add(new MockModel(diploid, new double[] {0.60, 0.25, 0.15}));
    models.add(new MockModel(diploid, new double[] {0.65, 0.15, 0.20}));
    final Family family = makeFamilySex(GenomeRelationships.SEX_MALE);
    final AbstractFamilyPosterior fp = getFamilyPosterior(priors, hdh, models, family);
    checkBest(fp.bestFather(), 1, -0.1286, 1.7263);
    checkBest(fp.bestMother(), 0, 0.1474, -0.1474);
    checkBest(fp.bestChild(0), 2, -0.3677, 0.9046);
    assertEquals(1, fp.numberChildren());
    assertTrue(fp.isInteresting());
    assertEquals(getExpectedNonIdentity("testDiploidDiploidDiploidDenovo"), fp.getNonIdentityPosterior(), 1e-4);
  }
  //all autosomes with denovo priors son should be called denovo
  // See testdiploiddiploiddiploid_denovo_2 in spreadsheet
  public void testDiploidDiploidDiploidDenovoSon() throws InvalidParamsException {
    final GenomePriorParams priors = new GenomePriorParamsBuilder().denovoRef(2.3E-3).denovoNonRef(2.3E-6).create();
//    final GenomePriorParams priors = new GenomePriorParamsBuilder().denovoRef(0).denovoNonRef(0).create();
    final PossibilityArithmetic arith = SimplePossibility.SINGLETON;
    final Description desc = new DescriptionCommon("A", "C");
    final HypothesesPrior<Description> none = new HypothesesNone<>(DescriptionNone.SINGLETON, arith, 0);
    final HypothesesPrior<Description> haploid = new MockHyp(desc, arith, true, new double[] {0.7, 0.3});
    final HypothesesPrior<Description> diploid = new MockHyp(desc, arith, false, new double[] {0.80, 0.10, 0.10});
    final HaploidDiploidHypotheses<HypothesesPrior<Description>> hdh = new HaploidDiploidHypotheses<>(none, haploid, diploid, false, null);
    final List<ModelInterface<?>> models = new ArrayList<>();
    models.add(new MockModel(diploid, new double[] {0.999, 0.0001, 0.0009}));
    models.add(new MockModel(diploid, new double[] {0.998, 0.00015, 0.00185}));
    models.add(new MockModel(diploid, new double[] {0.00005, 0.00005, 0.9999}));
    final Family family = makeFamilySex(GenomeRelationships.SEX_MALE);
    final AbstractFamilyPosterior fp = getFamilyPosterior(priors, hdh, models, family);
    checkBest(fp.bestFather(), 0, 2.9647, -2.9647);
    checkBest(fp.bestMother(), 0, 2.2431, -2.2431);
    checkBest(fp.bestChild(0), 2, 3.2981, 3.2981);
    if (!(this instanceof FamilyPosteriorFBTest)) { // XXX Come back here and implement denovo posterior for PPP
      assertEquals(1.5122, fp.bestChild(0).getDeNovoPosterior(), 1e-4);
    }
    assertEquals(1, fp.numberChildren());
    assertTrue(fp.isInteresting());
    assertEquals(getExpectedNonIdentity("testDiploidDiploidDiploidDenovoSon"), fp.getNonIdentityPosterior(), 1e-4);
  }
  //all autosomes with denovo priors son should be called denovo
  // See testdiploiddiploiddiploid_denovo_2 in spreadsheet
  public void testDiploidDiploidDiploidDenovoSonNoPrior() throws InvalidParamsException {
    final GenomePriorParams priors = new GenomePriorParamsBuilder().denovoRef(0).denovoNonRef(0).contraryProbability(1).create();
    final PossibilityArithmetic arith = SimplePossibility.SINGLETON;
    final Description desc = new DescriptionCommon("A", "C");
    final HypothesesPrior<Description> none = new HypothesesNone<>(DescriptionNone.SINGLETON, arith, 0);
    final HypothesesPrior<Description> haploid = new MockHyp(desc, arith, true, new double[] {0.7, 0.3});
    final HypothesesPrior<Description> diploid = new MockHyp(desc, arith, false, new double[] {0.80, 0.10, 0.10});
    final HaploidDiploidHypotheses<HypothesesPrior<Description>> hdh = new HaploidDiploidHypotheses<>(none, haploid, diploid, false, null);
    final List<ModelInterface<?>> models = new ArrayList<>();
    models.add(new MockModel(diploid, new double[] {0.999, 0.0001, 0.0009}));
    models.add(new MockModel(diploid, new double[] {0.998, 0.00015, 0.00185}));
    models.add(new MockModel(diploid, new double[] {0.00005, 0.00005, 0.9999}));
    final Family family = makeFamilySex(GenomeRelationships.SEX_MALE);
    final AbstractFamilyPosterior fp = getFamilyPosterior(priors, hdh, models, family);
    checkBest(fp.bestFather(), 0, 0.9867, -0.9867);
    checkBest(fp.bestMother(), 0, -0.1251, 0.1251);
    checkBest(fp.bestChild(0), 2, 1.4033, 1.4033);
    if (!(this instanceof FamilyPosteriorFBTest)) { // XXX Come back here and implement denovo posterior for PPP
      assertNull(fp.bestChild(0).getDeNovoPosterior());
    }
    assertEquals(1, fp.numberChildren());
    assertTrue(fp.isInteresting());
    assertEquals(getExpectedNonIdentity("testDiploidDiploidDiploidDenovoSonNoPrior"), fp.getNonIdentityPosterior(), 1e-4);
  }

  /**
   * @return a genome prior params with denovo priors zeroed
   * @throws InvalidParamsException
   */
  public GenomePriorParams getGenomePriorParams() throws InvalidParamsException {
    return new GenomePriorParamsBuilder().denovoRef(0.0).denovoNonRef(0.0).contraryProbability(1).create();
  }
}
