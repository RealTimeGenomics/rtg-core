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

package com.rtg.variant.bayes.multisample.family;


import static com.rtg.mode.DNARangeAT.A;
import static com.rtg.mode.DNARangeAT.C;
import static com.rtg.mode.DNARangeAT.G;
import static com.rtg.mode.DNARangeAT.T;

import java.util.ArrayList;
import java.util.List;

import com.rtg.relation.Family;
import com.rtg.util.InvalidParamsException;
import com.rtg.variant.GenomePriorParams;
import com.rtg.variant.GenomePriorParamsBuilder;
import com.rtg.variant.Variant;
import com.rtg.variant.VariantOutputLevel;
import com.rtg.variant.VariantParams;
import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.Model;
import com.rtg.variant.bayes.ModelInterface;
import com.rtg.variant.bayes.ModelNone;
import com.rtg.variant.bayes.ModelTest;
import com.rtg.variant.bayes.NoAlleleBalance;
import com.rtg.variant.bayes.multisample.HaploidDiploidHypotheses;
import com.rtg.variant.bayes.multisample.HypothesisScore;
import com.rtg.variant.bayes.snp.DescriptionSnp;
import com.rtg.variant.bayes.snp.EvidenceQ;
import com.rtg.variant.bayes.snp.HypothesesNone;
import com.rtg.variant.bayes.snp.HypothesesPrior;
import com.rtg.variant.bayes.snp.HypothesesSnp;
import com.rtg.variant.bayes.snp.StatisticsSnp;
import com.rtg.variant.util.arithmetic.SimplePossibility;

/**
 */
public class FamilyPosteriorTest extends AbstractFamilyPosteriorTest {

  @Override
  protected AbstractFamilyPosterior getFamilyPosterior(final GenomePriorParams priors, final HaploidDiploidHypotheses<HypothesesPrior<Description>> hdh, final List<ModelInterface<?>> models, final Family family) {
    for (ModelInterface<?> model : models) {
      model.freeze();
    }
    return new FamilyPosterior(family, priors, models, hdh);
  }

  protected GenomePriorParams mPriors = null;
  protected HypothesesPrior<Description> mHypotheses = null;

  @Override
  public void setUp() throws Exception {
    super.setUp();
    mPriors = getGenomePriorParams();
    mHypotheses = new HypothesesSnp(SimplePossibility.SINGLETON, mPriors, false, 0);
  }

  @Override
  protected void tearDown() throws Exception {
    super.tearDown();
    mPriors = null;
    mHypotheses = null;
  }

  protected AbstractFamilyPosterior getFamilyPosterior(List<ModelInterface<?>> models, Family family) {
    for (ModelInterface<?> model : models) {
      model.freeze();
    }
    return new FamilyPosterior(family, mPriors, models, new HaploidDiploidHypotheses<>(HypothesesNone.SINGLETON, null, mHypotheses));
  }

  ModelInterface<?> getModel() throws InvalidParamsException {
    return new Model<>(mHypotheses, new StatisticsSnp(mHypotheses.description()), new NoAlleleBalance());
  }

  /*
   * Father has A:A
   * Mother has A:A
   * Child has A:A
   */
  public void testEqual() throws InvalidParamsException {
    final Family family = FamilyCallerTest.makeFamily(FATHER, MOTHER, "c");
    final ModelInterface<?> father = getModel();
    final ModelInterface<?> mother = getModel();
    final ModelInterface<?> child = getModel();

    // 10 A reads for father, mother and child
    final EvidenceQ evidence = new EvidenceQ(DescriptionSnp.SINGLETON, 0, 0, 0, 0.1, 0.1, true, false, true, false, false);
    for (int i = 0; i < 10; ++i) {
      father.increment(evidence);
      mother.increment(evidence);
      child.increment(evidence);
    }
    final List<ModelInterface<?>> models = new ArrayList<>();
    models.add(father);
    models.add(mother);
    models.add(child);

    final AbstractFamilyPosterior fp = getFamilyPosterior(models, family);
    final HypothesisScore bestFather = fp.bestFather();
    final HypothesisScore bestMother = fp.bestMother();
    final HypothesisScore bestChild = fp.bestChild(0);
    assertEquals("A:A", mHypotheses.name(bestFather.hypothesis()));
    assertEquals("A:A", mHypotheses.name(bestMother.hypothesis()));
    assertEquals("A:A", mHypotheses.name(bestChild.hypothesis()));
    assertFalse(fp.isInteresting());

    /*
    System.err.println(fp.dumpMarginals());
    System.err.println(bestFather);
    System.err.println(bestMother);
    System.err.println(bestChild);
     */
  }

  /*
   * Father has A:A
   * Mother has A:C
   * Child has A:C
   */
  public void testCalculations1() throws InvalidParamsException {
    final Family family = FamilyCallerTest.makeFamily(FATHER, MOTHER, "c");
    final ModelInterface<?> father = getModel();
    final ModelInterface<?> mother = getModel();
    final ModelInterface<?> child = getModel();

    // 10 A reads for father
    for (int i = 0; i < 10; ++i) {
      father.increment(new EvidenceQ(DescriptionSnp.SINGLETON, 0, 0, 0, 0.1, 0.1, true, false, true, false, false));
    }

    // 5 A reads and 5 C reads for mother and child
    for (int i = 0; i < 5; ++i) {
      mother.increment(new EvidenceQ(DescriptionSnp.SINGLETON, 0, 0, 0, 0.1, 0.1, true, false, true, false, false));
      mother.increment(new EvidenceQ(DescriptionSnp.SINGLETON, 1, 0, 0, 0.1, 0.1, true, false, true, false, false));
      child.increment(new EvidenceQ(DescriptionSnp.SINGLETON, 0, 0, 0, 0.1, 0.1, true, false, true, false, false));
      child.increment(new EvidenceQ(DescriptionSnp.SINGLETON, 1, 0, 0, 0.1, 0.1, true, false, true, false, false));
    }

    final List<ModelInterface<?>> models = new ArrayList<>();
    models.add(father);
    models.add(mother);
    models.add(child);
    final AbstractFamilyPosterior fp = getFamilyPosterior(models, family);
    final HypothesisScore bestFather = fp.bestFather();
    final HypothesisScore bestMother = fp.bestMother();
    final HypothesisScore bestChild = fp.bestChild(0);
    assertEquals("A:A", mHypotheses.name(bestFather.hypothesis()));
    assertEquals("A:C", mHypotheses.name(bestMother.hypothesis()));
    assertEquals("A:C", mHypotheses.name(bestChild.hypothesis()));
    assertTrue(fp.isInteresting());

    /*
    System.err.println(fp.dumpMarginals());
    System.err.println(bestFather);
    System.err.println(bestMother);
    System.err.println(bestChild);
     */
  }

  AbstractFamilyPosterior makeFamily(String... members) {
    assert members.length > 2;
    final String[] children = new String[members.length - 2];
    System.arraycopy(members, 2, children, 0, children.length);
    final String[] childnames = new String[children.length];
    for (int i = 0; i < childnames.length; ++i) {
      childnames[i] = "c" + i;
    }
    final Family family = FamilyCallerTest.makeFamily(FATHER, MOTHER, childnames);

    final ArrayList<ModelInterface<?>> interfaces = new ArrayList<>();
    for (String member : members) {
      final ModelInterface<?> person = getModel();
      for (final char c : member.toCharArray()) {
        FamilyCallerTest.increment(person, "" + c);
      }
      interfaces.add(person);
    }
    return getFamilyPosterior(interfaces, family);
  }

  public void testNonIdentity() {
    final AbstractFamilyPosterior first = makeFamily("AAAAA", "AAAAA", "AAAAA");
    final AbstractFamilyPosterior second = makeFamily("AAAAA", "GGGGG", "AAGGG");

    assertTrue(first.getNonIdentityPosterior() < second.getNonIdentityPosterior());
    assertTrue(first.getNonIdentityPosterior() <= 0);
    assertTrue(second.getNonIdentityPosterior() >= 0);

  }

  public void testNonIdentityImproved() {
    final AbstractFamilyPosterior first = makeFamily("AAAAA", "GGGGG", "AAGGG");
    final AbstractFamilyPosterior second = makeFamily("AAAAA", "GGGGGGG", "AAGGG");
    assertTrue(first.getNonIdentityPosterior() < second.getNonIdentityPosterior());
  }
  public void testNonIdentityLessCertain() {
    final AbstractFamilyPosterior first = makeFamily("AAAAA", "AAGGG", "AAGGG");
    final AbstractFamilyPosterior second = makeFamily("AAAAA", "GGGGG", "AAGGG");
    assertTrue(first.getNonIdentityPosterior() < second.getNonIdentityPosterior());
  }
  public void testNonIdentityLessCertainChild() {
    final AbstractFamilyPosterior first = makeFamily("AAAAA", "AAGGG", "AAAAG");
    final AbstractFamilyPosterior second = makeFamily("AAAAA", "GGGGG", "AAGGG");
    assertTrue(first.getNonIdentityPosterior() < second.getNonIdentityPosterior());
  }
  public void testNonIdentityMultiChild() {
    final AbstractFamilyPosterior first = makeFamily("AAAAA", "AAGGG", "AAAAG", "AAAAA");
    final AbstractFamilyPosterior second = makeFamily("AAAAA", "GGGGG", "AAGGG", "AAAAAG");
    final AbstractFamilyPosterior third = makeFamily("AAAAA", "GGGGG", "AAGGGG", "AAAAA");
    assertTrue(first.getNonIdentityPosterior() < second.getNonIdentityPosterior());
    assertTrue(first.getNonIdentityPosterior() < third.getNonIdentityPosterior());
  }


  public void testIndividualNonIdentity() {
    final AbstractFamilyPosterior first = makeFamily("AAAAA", "AAAAA", "AAAAA");
    final AbstractFamilyPosterior second = makeFamily("AAAAA", "GGGGG", "AAGGG");
    final List<HypothesisScore> firstBest = bests(first);
    final List<HypothesisScore> secondBest = bests(second);
    for (int i = 0; i < firstBest.size(); ++i) {
      assertTrue(firstBest.get(i).nonIdentityPosterior() < secondBest.get(i).nonIdentityPosterior());
    }
  }

  public void testIndividualNonIdentityImproved() {
    final AbstractFamilyPosterior first = makeFamily("AAAAA", "GGGGG", "AAGGG");
    final AbstractFamilyPosterior second = makeFamily("AAAAA", "GGGGGG", "AAGGG");
    compareFamilyMembers(first, second, false, true, true);
  }

  public void testIndividualNonIdentityImprovedFather() {
    final AbstractFamilyPosterior first = makeFamily("AAAAA", "GGGGG", "AAGGG");
    final AbstractFamilyPosterior second = makeFamily("AAAAG", "GGGGG", "AAGGG");
    compareFamilyMembers(first, second, true, true, true);
  }

  public void testIndividualNonIdentityImprovedChild() {
    final AbstractFamilyPosterior first = makeFamily("AAAAA", "GGGGG", "AAAGG");
    final AbstractFamilyPosterior second = makeFamily("AAAAG", "GGGGG", "AAGGG");
    compareFamilyMembers(first, second, true, true, true);
  }

  public void testIndividualRandomness() {
    final AbstractFamilyPosterior first = makeFamily("ACGTACGT", "ACGTACGT", "ACGTACGT");
    final AbstractFamilyPosterior second = makeFamily("ACGTACGTACGT", "ACGTACGTACGT", "ACGTACGTACGT");
    compareFamilyMembers(first, second, true, true, true);
  }

  public void testIndividualMultiKids() {
    final AbstractFamilyPosterior first = makeFamily("AAAAA", "GGGGG", "AAGGG", "GGGGG");
    final AbstractFamilyPosterior second = makeFamily("AAAAA", "GGGGG", "AAGGG", "GGGGGAA");
    compareFamilyMembers(first, second, false, false, false, false);
  }
  public void testIndividualMultiKids2() {
    final AbstractFamilyPosterior first = makeFamily("AAAAA", "GGGGG", "AAGGG", "GGGGG");
    final AbstractFamilyPosterior second = makeFamily("AAAAA", "GGGGG", "AAGGG", "GGGGGG");
    compareFamilyMembers(first, second, true, true, true, true);
  }

  private void compareFamilyMembers(AbstractFamilyPosterior first, AbstractFamilyPosterior second, boolean... higher) {
    final List<HypothesisScore> firstBest = bests(first);
    final List<HypothesisScore> secondBest = bests(second);
    assertEquals(higher.length, firstBest.size());
    for (int i = 0; i < firstBest.size(); ++i) {
      final double firstScore = firstBest.get(i).nonIdentityPosterior();
      final double secondScore = secondBest.get(i).nonIdentityPosterior();
      if (higher[i]) {
        assertTrue("Individual " + i + " was not more confident: " + firstScore + " !< " + secondScore, firstScore < secondScore);
      } else {
        assertTrue("Individual " + i + " was more confident: " + firstScore + " !> " + secondScore, firstScore > secondScore);
      }
    }
  }

  List<HypothesisScore> bests(AbstractFamilyPosterior fp) {
    final ArrayList<HypothesisScore> b = new ArrayList<>();
    b.add(fp.bestFather());
    b.add(fp.bestMother());
    for (int i = 0; i < fp.mChildren.size(); ++i) {
      b.add(fp.bestChild(i));
    }
    return b;
  }
  /*
   * Father has A:A
   * Mother has A:C
   * Child has A:C
   * Child has A:A
   */
  public void testCalculations2() throws InvalidParamsException {
    final Family family = FamilyCallerTest.makeFamily(FATHER, MOTHER, "c", "c2");
    final ModelInterface<?> father = getModel();
    final ModelInterface<?> mother = getModel();
    final ModelInterface<?> child = getModel();
    final ModelInterface<?> child2 = getModel();

    // 10 A reads for father and child2
    for (int i = 0; i < 10; ++i) {
      father.increment(new EvidenceQ(DescriptionSnp.SINGLETON, 0, 0, 0, 0.1, 0.1, true, false, true, false, false));
      child2.increment(new EvidenceQ(DescriptionSnp.SINGLETON, 0, 0, 0, 0.1, 0.1, true, false, true, false, false));
    }

    // 5 A reads and 5 C reads for mother and child
    for (int i = 0; i < 5; ++i) {
      mother.increment(new EvidenceQ(DescriptionSnp.SINGLETON, 0, 0, 0, 0.1, 0.1, true, false, true, false, false));
      mother.increment(new EvidenceQ(DescriptionSnp.SINGLETON, 1, 0, 0, 0.1, 0.1, true, false, true, false, false));
      child.increment(new EvidenceQ(DescriptionSnp.SINGLETON, 0, 0, 0, 0.1, 0.1, true, false, true, false, false));
      child.increment(new EvidenceQ(DescriptionSnp.SINGLETON, 1, 0, 0, 0.1, 0.1, true, false, true, false, false));
    }

    final List<ModelInterface<?>> models = new ArrayList<>();
    models.add(father);
    models.add(mother);
    models.add(child);
    models.add(child2);
    final AbstractFamilyPosterior fp = getFamilyPosterior(models, family);
    final HypothesisScore bestFather = fp.bestFather();
    final HypothesisScore bestMother = fp.bestMother();
    final HypothesisScore bestChild1 = fp.bestChild(0);
    final HypothesisScore bestChild2 = fp.bestChild(1);
    //    System.err.println(fp.dumpMarginals());
    //    System.err.println(bestFather);
    //    System.err.println(bestMother);
    //    System.err.println(bestChild1);
    //    System.err.println(bestChild2);
    assertEquals("A:A", mHypotheses.name(bestFather.hypothesis()));
    assertEquals("A:C", mHypotheses.name(bestMother.hypothesis()));
    assertEquals("A:C", mHypotheses.name(bestChild1.hypothesis()));
    assertEquals("A:A", mHypotheses.name(bestChild2.hypothesis()));
    assertTrue(fp.isInteresting());
  }
  /*
   * Father has A:C
   * Mother has C:G
   * Child has A:C
   * Child has C:G
   */
  public void testFullCats() throws InvalidParamsException {
    final Family family = FamilyCallerTest.makeFamily(FATHER, MOTHER, "c", "c2");
    final ModelInterface<?> father = getModel();
    final ModelInterface<?> mother = getModel();
    final ModelInterface<?> child = getModel();
    final ModelInterface<?> child2 = getModel();
    // 10 A reads and 10 C reads for father and child1
    // 10 C reads and 10 G reads for mother and child2
    for (int i = 0; i < 5; ++i) {
      father.increment(new EvidenceQ(DescriptionSnp.SINGLETON, A, 0, 0, 0.1, 0.1, true, false, true, false, false));
      father.increment(new EvidenceQ(DescriptionSnp.SINGLETON, C, 0, 0, 0.1, 0.1, true, false, true, false, false));

      child.increment(new EvidenceQ(DescriptionSnp.SINGLETON, A, 0, 0, 0.1, 0.1, true, false, true, false, false));
      child.increment(new EvidenceQ(DescriptionSnp.SINGLETON, C, 0, 0, 0.1, 0.1, true, false, true, false, false));

      mother.increment(new EvidenceQ(DescriptionSnp.SINGLETON, C, 0, 0, 0.1, 0.1, true, false, true, false, false));
      mother.increment(new EvidenceQ(DescriptionSnp.SINGLETON, G, 0, 0, 0.1, 0.1, true, false, true, false, false));

      child2.increment(new EvidenceQ(DescriptionSnp.SINGLETON, C, 0, 0, 0.1, 0.1, true, false, true, false, false));
      child2.increment(new EvidenceQ(DescriptionSnp.SINGLETON, G, 0, 0, 0.1, 0.1, true, false, true, false, false));
    }
    // Add some noise
    father.increment(new EvidenceQ(DescriptionSnp.SINGLETON, G, 0, 0, 0.1, 0.1, true, false, true, false, false));
    mother.increment(new EvidenceQ(DescriptionSnp.SINGLETON, T, 0, 0, 0.1, 0.1, true, false, true, false, false));
    child.increment(new EvidenceQ(DescriptionSnp.SINGLETON, T, 0, 0, 0.1, 0.1, true, false, true, false, false));
    child2.increment(new EvidenceQ(DescriptionSnp.SINGLETON, A, 0, 0, 0.1, 0.1, true, false, true, false, false));

    final List<ModelInterface<?>> models = new ArrayList<>();
    models.add(father);
    models.add(mother);
    models.add(child);
    models.add(child2);
    final AbstractFamilyPosterior fp = getFamilyPosterior(models, family);
    final HypothesisScore bestFather = fp.bestFather();
    final HypothesisScore bestMother = fp.bestMother();
    final HypothesisScore bestChild1 = fp.bestChild(0);
    final HypothesisScore bestChild2 = fp.bestChild(1);
    //    System.err.println("========");
    //    System.err.println("AllCats");
    //    System.err.println(fp.dumpMarginals());
    //    System.err.println(bestFather);
    //    System.err.println(bestMother);
    //    System.err.println(bestChild1);
    //    System.err.println(bestChild2);
    assertEquals("A:C", mHypotheses.name(bestFather.hypothesis()));
    assertEquals("C:G", mHypotheses.name(bestMother.hypothesis()));
    assertEquals("A:C", mHypotheses.name(bestChild1.hypothesis()));
    assertEquals("C:G", mHypotheses.name(bestChild2.hypothesis()));
    assertTrue(fp.isInteresting());
  }

  private ModelInterface<?> createModel(int aCount, int cCount) {
    final ModelInterface<?> childModel = getModel();
    for (int i = 0; i < aCount; ++i) {
      childModel.increment(new EvidenceQ(DescriptionSnp.SINGLETON, A, 0, 0, 0.1, 0.2, true, false, true, false, false));
    }

    for (int i = 0; i < cCount; ++i) {
      childModel.increment(new EvidenceQ(DescriptionSnp.SINGLETON, C, 0, 0, 0.1, 0.2, true, false, true, false, false));
    }
    childModel.freeze();
    return childModel;
  }

  public void testBorderChildCase() throws InvalidParamsException {
    final int readCount = 30;
    final int callThreshold = getCallThreshold(readCount);
    // callThreshold is now the number of C reads required to convert an A call to A:C


    final ModelInterface<?> father = getModel();
    final ModelInterface<?> mother = getModel();
    final ModelInterface<?> child = getModel();
    for (int i = 0; i < readCount; ++i) {
      father.increment(new EvidenceQ(DescriptionSnp.SINGLETON, A, 0, 0, 0.1, 0.2, true, false, true, false, false));

      mother.increment(new EvidenceQ(DescriptionSnp.SINGLETON, C, 0, 0, 0.1, 0.2, true, false, true, false, false));

      child.increment(new EvidenceQ(DescriptionSnp.SINGLETON, A, 0, 0, 0.1, 0.2, true, false, true, false, false));
    }
    //incrementCats(mother, new EvidenceQ(DescriptionSnp.SINGLETON, 0.1, A, 0.2));


    // with family calling should be able to do this with fewer C reads
    for (int i = 0; i < callThreshold / 5; ++i) {
      child.increment(new EvidenceQ(DescriptionSnp.SINGLETON, C, 0, 0, 0.1, 0.2, true, false, true, false, false));
    }

    final List<ModelInterface<?>> models = new ArrayList<>();
    models.add(father);
    models.add(mother);
    models.add(child);
    final Family family = FamilyCallerTest.makeFamily(FATHER, MOTHER, "c");
    final AbstractFamilyPosterior fp = getFamilyPosterior(models, family);
    final HypothesisScore bestFather = fp.bestFather();
    final HypothesisScore bestMother = fp.bestMother();
    final HypothesisScore bestChild1 = fp.bestChild(0);
    //      System.err.println("========");
    //      System.err.println("ChildBorder");
    //      System.err.println("Snp caller: " + bc);
    //      System.err.println(" Max reads: " + (callThreshold - 1));
    //
    //      System.err.println(fp.dumpMarginals());
    //      System.err.println(bestFather);
    //      System.err.println(bestMother);
    //      System.err.println(bestChild1);
    assertEquals("A:A", mHypotheses.name(bestFather.hypothesis()));
    assertEquals("C:C", mHypotheses.name(bestMother.hypothesis()));
    assertEquals("A:C", mHypotheses.name(bestChild1.hypothesis()));
  }

  public void testBorderParentCase() throws InvalidParamsException {
    final int readCount = 50;
    final int callThreshold = getCallThreshold(readCount);
    // callThreshold is now the number of C reads required to convert an A call to A:C

    final ModelInterface<?> father = getModel();
    final ModelInterface<?> mother = getModel();
    final ModelInterface<?> child = getModel();
    for (int i = 0; i < readCount; ++i) {
      father.increment(new EvidenceQ(DescriptionSnp.SINGLETON, A, 0, 0, 0.1, 0.2, true, false, true, false, false));

      mother.increment(new EvidenceQ(DescriptionSnp.SINGLETON, A, 0, 0, 0.1, 0.2, true, false, true, false, false));

      child.increment(new EvidenceQ(DescriptionSnp.SINGLETON, i % 2 == 0 ? A : C, 0, 0, 0.1, 0.2, true, false, true, false, false));
    }


    // with family calling should be able to do this with fewer C reads
    for (int i = 0; i < callThreshold; ++i) {
      mother.increment(new EvidenceQ(DescriptionSnp.SINGLETON, C, 0, 0, 0.1, 0.2, true, false, true, false, false));
    }

    final List<ModelInterface<?>> models = new ArrayList<>();
    models.add(father);
    models.add(mother);
    models.add(child);
    final Family family = FamilyCallerTest.makeFamily(FATHER, MOTHER, "c");
    final AbstractFamilyPosterior fp = getFamilyPosterior(models, family);
    final HypothesisScore bestFather = fp.bestFather();
    final HypothesisScore bestMother = fp.bestMother();
    final HypothesisScore bestChild1 = fp.bestChild(0);
    //      System.err.println("========");
    //      System.err.println("ParentBorder");
    //      System.err.println("Snp caller: " + bc);
    //      System.err.println(" Max reads: " + (callThreshold - 1));
    //
    //      System.err.println(fp.dumpMarginals());
    //      System.err.println(bestFather);
    //      System.err.println(bestMother);
    //      System.err.println(bestChild1);
    assertEquals("A:A", mHypotheses.name(bestFather.hypothesis()));
    assertEquals("A:C", mHypotheses.name(bestMother.hypothesis()));
    assertEquals("A:C", mHypotheses.name(bestChild1.hypothesis()));
  }

  private int getCallThreshold(int readCount) {
    int callThreshold = 0;
    while (true) {
      final ModelInterface<?> childModel = createModel(readCount, callThreshold);
      final Variant v = ModelTest.makeCalls(childModel, "blah", 0, 1, new byte[] {}, VariantParams.builder().callLevel(VariantOutputLevel.ALL).create());
      if (!v.getSample(0).getName().equals("A:A")) {
        break;
      }
      ++callThreshold;
    }
    return callThreshold;
  }

  public void testDenovoPosterior() {
    mPriors = new GenomePriorParamsBuilder().contraryProbability(1).create();
    final AbstractFamilyPosterior first = makeFamily("AAAAAAAAAAAAAAAAAA", "AAAAAAAAAAAAAAAAAAA", "AAAAAAAGGGGGGGGGGGGGGGGGGGGGGGGG");
    final AbstractFamilyPosterior second = makeFamily("AAAAAAAAAAAAAAAAAA", "AAAAAAAAAAAAAAAAAAAG", "AAAAAAAAAAGGGGGGGGGGGGGGGGGGGGGGGGG");
    final HypothesisScore firstScore = first.bestChild(0);
    final HypothesisScore secondScore = second.bestChild(0);
    assertTrue(firstScore.getDeNovoPosterior() > secondScore.getDeNovoPosterior());
  }

  public void testPolyploid() {

  }

  public void testAllNone() {
    final ModelInterface<?> father = ModelNone.SINGLETON;
    final ModelInterface<?> mother = ModelNone.SINGLETON;
    final ModelInterface<?> child = ModelNone.SINGLETON;
    final List<ModelInterface<?>> models = new ArrayList<>();
    models.add(father);
    models.add(mother);
    models.add(child);
    final Family family = FamilyCallerTest.makeFamily(FATHER, MOTHER, "c");
    final AbstractFamilyPosterior fp = getFamilyPosterior(models, family);
    assertNull(fp.bestFather());
    assertNull(fp.bestMother());
    assertNull(fp.bestChild(0));
  }
}

