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

import java.io.BufferedReader;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.List;

import com.rtg.relation.Family;
import com.rtg.relation.RelationshipsFileParser;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.StringUtils;
import com.rtg.variant.GenomePriorParams;
import com.rtg.variant.GenomePriorParamsBuilder;
import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.Model;
import com.rtg.variant.bayes.ModelInterface;
import com.rtg.variant.bayes.StatisticsInt;
import com.rtg.variant.bayes.snp.DescriptionSnp;
import com.rtg.variant.bayes.snp.HypothesesMock;
import com.rtg.variant.util.arithmetic.SimplePossibility;

import junit.framework.TestCase;

/**
 */
public class DiseasedFamilyPosteriorTest extends TestCase {
  protected double diseasePrior() {
    return 0.999;
  }

  protected DiseasedFamilyPosterior getDiseasedFamily(Family family, final int refNt, final List<ModelInterface<?>> models) {
    try {
      final GenomePriorParams priors = new GenomePriorParamsBuilder().create();
      final HypothesesDisease diseaseHypotheses = new HypothesesDisease(models.get(0).description(), diseasePrior(), refNt);
      final DiseasedFamilyPosterior diseasedFamilyPosterior = new DiseasedFamilyPosterior(priors, family, models.get(0).hypotheses(), diseaseHypotheses, models);
      diseasedFamilyPosterior.process();
      return diseasedFamilyPosterior;
    } catch (final InvalidParamsException e) {
      throw new RuntimeException(e);
    }
  }

  /** mother and second child are diseased. */
  public static final String RELATIONS = ""
    + "#parents" + StringUtils.LS
    + "#twins" + StringUtils.LS
    + "#cancer" + StringUtils.LS
    + "" + StringUtils.LS
    + "genome mother disease=" + true + " sex=female" + StringUtils.LS
    + "genome father disease=" + false + " sex=male" + StringUtils.LS
    + "genome twina\t\tdisease=" + false + " sex=female" + StringUtils.LS
    + "genome twinb\t\tdisease=" + true + StringUtils.LS
    + "parent-child\tfather\ttwina" + StringUtils.LS
    + "parent-child\t father\ttwinb" + StringUtils.LS
    + "parent-child  mother twina" + StringUtils.LS
    + "parent-child mother twinb" + StringUtils.LS
    + "original-derived father fathercancer contamination=0.03" + StringUtils.LS;

  /** Only mother diseased. */
  public static final String RELATIONS1 = ""
      + "genome father sex=male disease=" + false + StringUtils.LS
      + "genome mother sex=female disease=" + true + StringUtils.LS
      + "genome child disease=" + false + StringUtils.LS
      + "parent-child father child" + StringUtils.LS
      + "parent-child  mother child" + StringUtils.LS;

  /** Mother and child diseased. */
  public static final String RELATIONS2 = ""
      + "genome father sex=male disease=" + false + StringUtils.LS
      + "genome mother disease=" + true + StringUtils.LS
      + "genome child disease=" + true + StringUtils.LS
      + "parent-child father child" + StringUtils.LS
      + "parent-child  mother child" + StringUtils.LS;

  /** Mother and second child diseased. */
  public static final String RELATIONS3 = ""
      + "genome father sex=male disease=" + false + StringUtils.LS
      + "genome mother sex=female disease=" + true + StringUtils.LS
      + "genome childa disease=" + false + StringUtils.LS
      + "genome childb disease=" + true + StringUtils.LS
      + "parent-child father childa" + StringUtils.LS
      + "parent-child  mother childa" + StringUtils.LS
      + "parent-child father childb" + StringUtils.LS
      + "parent-child  mother childb" + StringUtils.LS;

  public void test1() throws Exception {
    final GenomePriorParams priors = new GenomePriorParamsBuilder().create();
    final Family f = Family.getFamily(RelationshipsFileParser.load(new BufferedReader(new StringReader(RELATIONS1))));
    assertEquals("mother", f.getMother());
    assertEquals("father", f.getFather());
    final List<ModelInterface<?>> models =  FamilyCallerTest.buildFamily(priors, 1, "AAAAAAAAAAA", "AAAAAAATTTTTTTT", "AAAAAAAATTTTTTTTTTT");
    final ModelInterface<?> model = models.get(0);
    final DiseasedFamilyPosterior diseased = getDiseasedFamily(f, 0, models);
    //System.err.println(dumpMarginals(diseased));

    assertEquals(0, diseased.bestDisease().hypothesis());
    assertEquals("A:T", model.name(diseased.bestMother().hypothesis()));
    assertEquals("A:A", model.name(diseased.bestFather().hypothesis()));
    assertEquals("A:T", model.name(diseased.bestChild(0).hypothesis()));
    assertFalse(diseased.isInteresting());
  }

  public void test3() throws Exception {
    final GenomePriorParams priors = new GenomePriorParamsBuilder().create();
    final Family f = Family.getFamily(RelationshipsFileParser.load(new BufferedReader(new StringReader(RELATIONS3))));
    assertEquals("mother", f.getMother());
    assertEquals("father", f.getFather());
    final List<ModelInterface<?>> models =  FamilyCallerTest.buildFamily(priors, 1, "AAAAAAAAAAA", "AAAAAAATTTTTTTT", "AAAAAAAAAAAAAAA", "AAAAAAAATTTTTTTTTTT");
    final ModelInterface<?> model = models.get(0);
    final DiseasedFamilyPosterior diseased = getDiseasedFamily(f, 0, models);
    //System.err.println(dumpMarginals(diseased));

    assertEquals(4, diseased.bestDisease().hypothesis());
    assertEquals("A:T", model.name(diseased.bestMother().hypothesis()));
    assertEquals("A:A", model.name(diseased.bestFather().hypothesis()));
    assertEquals("A:A", model.name(diseased.bestChild(0).hypothesis()));
    assertEquals("A:T", model.name(diseased.bestChild(1).hypothesis()));
    assertTrue(diseased.isInteresting());
  }

  // mother and second child are diseased
  public void test() throws Exception {
    final GenomePriorParams priors = new GenomePriorParamsBuilder().create();
    final Family f = Family.getFamily(RelationshipsFileParser.load(new BufferedReader(new StringReader(RELATIONS))));
    assertEquals("mother", f.getMother());
    assertEquals("father", f.getFather());
    final List<ModelInterface<?>> models =  FamilyCallerTest.buildFamily(priors, 1, "AAAAAAAAAAA", "AAAAAAATTTTTTTT", "AAAAAAAAAAAAAAA", "AAAAAAAATTTTTTTTTTT");
    final ModelInterface<?> model = models.get(0);
    final DiseasedFamilyPosterior diseased = getDiseasedFamily(f, 0, models);

    assertEquals(4, diseased.bestDisease().hypothesis());
    assertEquals("A:T", model.name(diseased.bestMother().hypothesis()));
    assertEquals("A:A", model.name(diseased.bestFather().hypothesis()));
    assertEquals("A:A", model.name(diseased.bestChild(0).hypothesis()));
    assertEquals("A:T", model.name(diseased.bestChild(1).hypothesis()));

    assertTrue(diseased.isInteresting());
  }

  public void testNoDiseased() throws Exception {
    final GenomePriorParams priors = new GenomePriorParamsBuilder().create();
    final Family f = Family.getFamily(RelationshipsFileParser.load(new BufferedReader(new StringReader(RELATIONS))));
    assertEquals("mother", f.getMother());
    assertEquals("father", f.getFather());
    assertFalse(f.isDiseased(f.getFather()));
    assertTrue(f.isDiseased(f.getMother()));
    assertFalse(f.isDiseased(f.getChildren()[0]));
    assertTrue(f.isDiseased(f.getChildren()[1]));
    final List<ModelInterface<?>> models =  FamilyCallerTest.buildFamily(priors, 1
            , "AAAAAAAAAAAAAA"
            , "AAAAAAAATTTTTTTT"
            , "AAAAAAAATTTTTTTT"
            , "AAAAAAAATTTTTTTT"
    );
    final ModelInterface<?> model = models.get(0);
    final DiseasedFamilyPosterior diseased = getDiseasedFamily(f, 0, models);
    //System.err.println(dumpMarginals(diseased));

    assertEquals(0, diseased.bestDisease().hypothesis());

    assertEquals(13.0569, diseased.bestDisease().genotypeMeasure().bestPosterior(), 1E-4);

    assertEquals("A:A", model.name(diseased.bestFather().hypothesis()));
    assertEquals("A:T", model.name(diseased.bestMother().hypothesis()));
    assertEquals("A:T", model.name(diseased.bestChild(0).hypothesis()));
    assertEquals("A:T", model.name(diseased.bestChild(1).hypothesis()));
    assertFalse(diseased.isInteresting());
  }

  // mother and first child are diseased
  public void testOneChild() throws Exception {
    final GenomePriorParams priors = new GenomePriorParamsBuilder().create();
    final Family f = Family.getFamily(RelationshipsFileParser.load(new BufferedReader(new StringReader(RELATIONS2))));
    final List<ModelInterface<?>> models =  FamilyCallerTest.buildFamily(priors, 1, "AAAAAAAAAAA", "AAAAAAATTTTTTTT", "AAAAAAAATTTTTTTTTTT");
    final ModelInterface<?> model = models.get(0);
    final DiseasedFamilyPosterior diseased = getDiseasedFamily(f, 0, models);
    //System.err.println(dumpMarginals(diseased));

    assertEquals(4, diseased.bestDisease().hypothesis());
    assertEquals("A:T", model.name(diseased.bestMother().hypothesis()));
    assertEquals("A:A", model.name(diseased.bestFather().hypothesis()));
    assertEquals("A:T", model.name(diseased.bestChild(0).hypothesis()));
    assertTrue(diseased.isInteresting());
  }

  private static class MockHyp extends HypothesesMock<Description> {
    MockHyp(final int refHyp) {
      super(DescriptionSnp.SINGLETON, SimplePossibility.SINGLETON, false, refHyp);
    }
    @Override
    protected void initPriors(double[] priors) {
      for (int i = 0; i < priors.length; i++) {
        priors[i] = arithmetic().one();
      }
    }
  }

  private static final MockHyp HYPOTHESES0 = new MockHyp(0);
  private static final MockHyp HYPOTHESES1 = new MockHyp(1);

  private static class MockModel extends Model<Description> {
    MockModel(final double[] posteriors) {
      super(HYPOTHESES0, new StatisticsInt(HYPOTHESES0.description()));
      assert posteriors.length == mPosteriors.length;
      for (int i = 0; i < posteriors.length; i++) {
        mPosteriors[i] = arithmetic().prob2Poss(posteriors[i]);
      }
    }
  }

  private static class MockModelC extends Model<Description> {
    MockModelC(final double[] posteriors) {
      super(HYPOTHESES1, new StatisticsInt(HYPOTHESES1.description()));
      assert posteriors.length == mPosteriors.length;
      for (int i = 0; i < posteriors.length; i++) {
        mPosteriors[i] = arithmetic().prob2Poss(posteriors[i]);
      }
    }
  }

  private ModelInterface<?> make(final double[] prior) {
    return new MockModel(prior);
  }

  private ModelInterface<?> makeC(final double[] prior) {
    return new MockModelC(prior);
  }

  private static final double[] POST1 = new double[10];
  static {
    POST1[0] = 0.5; //A
    POST1[HYPOTHESES0.code().code(0, 3)] = 0.5; //A:T
  }

  private static final double[] POST2 = new double[10];
  static {
    POST2[3] = 0.5; //T
    POST2[HYPOTHESES0.code().code(0, 3)] = 0.5; //A:T
  }

  private static final double[] POST1C = new double[10];
  static {
    POST1C[0] = 0.5; //A
    POST1C[HYPOTHESES1.code().code(0, 3)] = 0.5; //A:T
  }

  private static final double[] POST2C = new double[10];
  static {
    POST2C[3] = 0.5; //T
    POST2C[HYPOTHESES1.code().code(0, 3)] = 0.5; //A:T
  }

  //
  public void testOneDiseasedChildReduced() throws Exception {
    //see DiseasedFamilyTest.xslx for details of numbers
    final Family family = Family.getFamily(RelationshipsFileParser.load(new BufferedReader(new StringReader(RELATIONS2))));
    final ModelInterface<?> father = make(POST1);
    final ModelInterface<?> mother = make(POST2);
    final ModelInterface<?> child = make(POST1);
    final List<ModelInterface<?>> list = new ArrayList<>();
    list.add(father);
    list.add(mother);
    list.add(child);
    final DiseasedFamilyPosterior diseased = getDiseasedFamily(family, 0, list);
    //System.err.println(dumpMarginals(diseased));
    assertEquals(0, diseased.bestDisease().hypothesis());
    assertEquals(2.281874459, diseased.bestDisease().genotypeMeasure().bestPosterior(), 0.00001);

    assertEquals("A:T", father.name(diseased.bestMother().hypothesis()));
    assertEquals(1.386294, diseased.bestMother().genotypeMeasure().bestPosterior(), 0.0001);

    assertEquals("A:T", father.name(diseased.bestFather().hypothesis()));
    assertEquals(1.386294, diseased.bestFather().genotypeMeasure().bestPosterior(), 0.0001);

    assertEquals("A:T", father.name(diseased.bestChild(0).hypothesis()));
    assertEquals(0.405455, diseased.bestChild(0).genotypeMeasure().bestPosterior(), 0.0001);

    assertFalse(diseased.isInteresting());
  }

  //
  public void testOneHealthyChildReduced() throws Exception {
    //see DiseasedFamilyTest.xslx for details of numbers
    final Family family = Family.getFamily(RelationshipsFileParser.load(new BufferedReader(new StringReader(RELATIONS1))));
    final ModelInterface<?> father = make(POST1);
    final ModelInterface<?> mother = make(POST2);
    final ModelInterface<?> child = make(POST1);
    final List<ModelInterface<?>> list = new ArrayList<>();
    list.add(father);
    list.add(mother);
    list.add(child);
    final DiseasedFamilyPosterior diseased = getDiseasedFamily(family, 0, list);
    //System.err.println(dumpMarginals(diseased));
    assertEquals(0, diseased.bestDisease().hypothesis());
    assertEquals(3.157343197, diseased.bestDisease().genotypeMeasure().bestPosterior(), 0.00002);

    assertEquals("A:T", father.name(diseased.bestMother().hypothesis()));
    assertEquals(0.6931, diseased.bestMother().genotypeMeasure().bestPosterior(), 0.0001);

    assertEquals("A:T", father.name(diseased.bestFather().hypothesis()));
    assertEquals(0.6931, diseased.bestFather().genotypeMeasure().bestPosterior(), 0.0001);

    assertEquals("A:T", father.name(diseased.bestChild(0).hypothesis()));
    assertEquals(1.6094, diseased.bestChild(0).genotypeMeasure().bestPosterior(), 0.0001);

    assertFalse(diseased.isInteresting());
  }

  //swap around nucleotides A->C T->G
  public void testOneHealthyChildReducedX() throws Exception {
    //regression test between fast and slow
    final Family family = Family.getFamily(RelationshipsFileParser.load(new BufferedReader(new StringReader(RELATIONS1))));
    final ModelInterface<?> father = makeC(POST1C);
    final ModelInterface<?> mother = makeC(POST2C);
    final ModelInterface<?> child = makeC(POST1C);
    final List<ModelInterface<?>> list = new ArrayList<>();
    list.add(father);
    list.add(mother);
    list.add(child);
    final DiseasedFamilyPosterior diseased = getDiseasedFamily(family, 1 /*refNt=C*/, list);
    //System.err.println(dumpMarginals(diseased));

    assertEquals(0, diseased.bestDisease().hypothesis());
    assertEquals(2.48345056530216, diseased.bestDisease().genotypeMeasure().bestPosterior(), 0.0001);

    assertEquals("A:T", father.name(diseased.bestMother().hypothesis()));
    assertEquals(0.6931, diseased.bestMother().genotypeMeasure().bestPosterior(), 0.0001);

    assertEquals("A:T", father.name(diseased.bestFather().hypothesis()));
    assertEquals(0.6931, diseased.bestFather().genotypeMeasure().bestPosterior(), 0.0001);

    assertEquals("A:T", father.name(diseased.bestChild(0).hypothesis()));
    assertEquals(1.6094, diseased.bestChild(0).genotypeMeasure().bestPosterior(), 0.0001);

    assertFalse(diseased.isInteresting());
  }

  public void testNonMendelian() throws Exception {
    final GenomePriorParams priors = new GenomePriorParamsBuilder().create();
    final Family f = Family.getFamily(RelationshipsFileParser.load(new BufferedReader(new StringReader(DiseasedFamilyPosteriorTest.RELATIONS))));
    assertEquals("father", f.getFather());
    assertEquals("mother", f.getMother());
    final List<ModelInterface<?>> models =  FamilyCallerTest.buildFamily(priors, 1
            , "AAAAAAAAAAAAAA"
            , "AAAAAAAATTTTTTTT"
            , "AAAAAAAATTTTTTTT"
            , "GGGGGGGGGGGGGGGG"
    );
    final ModelInterface<?> model = models.get(0);
    final DiseasedFamilyPosterior diseased = getDiseasedFamily(f, 0, models);
    //System.err.println(dumpMarginals(diseased));
    assertEquals(0, diseased.bestDisease().hypothesis());

    assertEquals("A:G", model.name(diseased.bestFather().hypothesis()));
    assertEquals("A:T", model.name(diseased.bestMother().hypothesis()));
    assertEquals("A:T", model.name(diseased.bestChild(0).hypothesis()));
    assertEquals("G:T", model.name(diseased.bestChild(1).hypothesis()));
    assertFalse(diseased.isInteresting());
  }

}
