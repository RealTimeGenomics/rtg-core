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


import java.util.ArrayList;
import java.util.List;

import com.rtg.relation.Family;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.PortableRandom;
import com.rtg.variant.GenomePriorParams;
import com.rtg.variant.GenomePriorParamsBuilder;
import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.ModelInterface;
import com.rtg.variant.bayes.multisample.HaploidDiploidHypotheses;
import com.rtg.variant.bayes.multisample.HypothesisScore;
import com.rtg.variant.bayes.snp.DescriptionSnp;
import com.rtg.variant.bayes.snp.EvidenceQ;
import com.rtg.variant.bayes.snp.HypothesesNone;
import com.rtg.variant.bayes.snp.HypothesesPrior;

/**
 */
public class FastFamilyPosteriorTest extends FamilyPosteriorTest {

  @Override
  protected AbstractFamilyPosterior getFamilyPosterior(List<ModelInterface<?>> models, Family family) {
    final GenomePriorParams priors = mPriors;
    final List<ModelInterface<?>> list = new ArrayList<>(models);
    for (ModelInterface<?> model : models) {
      model.freeze();
    }
    return new FastFamilyPosterior(family, priors, list, new HaploidDiploidHypotheses<>(HypothesesNone.SINGLETON, null, mHypotheses));
  }

  @Override
  protected AbstractFamilyPosterior getFamilyPosterior(final GenomePriorParams priors, final HaploidDiploidHypotheses<HypothesesPrior<Description>> hdh, final List<ModelInterface<?>> models, final Family family) {
    for (ModelInterface<?> model : models) {
      model.freeze();
    }
    return new FastFamilyPosterior(family, priors, models, hdh);
  }

  private void checkRandom(final PortableRandom r) throws InvalidParamsException {
    final ModelInterface<?> father = getModel();
    final ModelInterface<?> mother = getModel();
    final List<ModelInterface<?>> models = new ArrayList<>();
    final int numChildren = 1 + r.nextInt(5);
    models.add(father);
    models.add(mother);
    final List<String> children = new ArrayList<>();
    for (int c = 0; c < numChildren; ++c) {
      models.add(getModel());
      children.add("c" + c);
    }
    // Choose genotypes for mother and father
    final int fa = r.nextInt(4);
    final int fb = r.nextInt(4);
    final int ma = r.nextInt(4);
    final int mb = r.nextInt(4);
    // Choose read counts for mother and father
    final int fc = 1 + r.nextInt(20);
    final int mc = 1 + r.nextInt(20);
    // Add reads for father
    for (int i = 0; i < fc; ++i) {
      //incrementCats(father, new ProbabilityQ(r.nextBoolean() ? fa : fb, r.nextDouble()));
      father.increment(new EvidenceQ(DescriptionSnp.SINGLETON, r.nextBoolean() ? fa : fb, 0, 0, 0.1, r.nextDouble(), true, false, true, false, false));
    }
    father.freeze();
    // Add reads for mother
    for (int i = 0; i < mc; ++i) {
      //incrementCats(mother, new ProbabilityQ(r.nextBoolean() ? ma : mb, r.nextDouble()));
      mother.increment(new EvidenceQ(DescriptionSnp.SINGLETON, r.nextBoolean() ? ma : mb, 0, 0, 0.1, r.nextDouble(), true, false, true, false, false));
    }
    mother.freeze();
    // Note: this always does Mendelian valid children
    for (int i = 0; i < numChildren; ++i) {
      // Choose genotype of child
      final int ca = r.nextBoolean() ? fa : fb;
      final int cb = r.nextBoolean() ? ma : mb;
      final int cc = 1 + r.nextInt(20);
      // Add reads for child
      final ModelInterface<?> child = models.get(i + Family.FIRST_CHILD_INDEX);
      for (int j = 0; j < cc; ++j) {
        //incrementCats(child, new ProbabilityQ(r.nextBoolean() ? ca : cb, r.nextDouble()));
        child.increment(new EvidenceQ(DescriptionSnp.SINGLETON, r.nextBoolean() ? ca : cb, 0, 0, 0.1, r.nextDouble(), true, false, true, false, false));
      }
      child.freeze();
    }

    final GenomePriorParams priors = getGenomePriorParams();
    final Family family = FamilyCallerTest.makeFamily("f", "m", children.toArray(new String[0]));
    final AbstractFamilyPosterior slow = new FamilyPosterior(family, priors, models, new HaploidDiploidHypotheses<>(HypothesesNone.SINGLETON, null, mHypotheses));
    final AbstractFamilyPosterior fast = new FastFamilyPosterior(family, priors, models, new HaploidDiploidHypotheses<>(HypothesesNone.SINGLETON, null, mHypotheses));
    assertEquals(fast.bestFather().hypothesis(), slow.bestFather().hypothesis());
    assertEquals(fast.bestFather().genotypeMeasure().bestPosterior(), slow.bestFather().genotypeMeasure().bestPosterior(), 1E-2);

    assertEquals(fast.bestMother().hypothesis(), slow.bestMother().hypothesis());
    assertEquals(fast.bestMother().genotypeMeasure().bestPosterior(), slow.bestMother().genotypeMeasure().bestPosterior(), 1E-2);

    for (int i = 0; i < numChildren; ++i) {
      final HypothesisScore f = fast.bestChild(i);
      final HypothesisScore s = slow.bestChild(i);
      assertEquals(f.hypothesis(), s.hypothesis());
      assertEquals(f.genotypeMeasure().bestPosterior(), s.genotypeMeasure().bestPosterior(), 1E-3);
    }
  }

  public void testBradyBunch() throws Exception {
    final PortableRandom r = new PortableRandom(42);
    for (int k = 0; k < 10; ++k) {
      checkRandom(r);
    }
  }
  public void testFamilyScalingNonDenovo() throws Exception {
    // test case inspired by calls where all children are "de novo"
    // it isn't clear what the answer should be here.
    mPriors = new GenomePriorParamsBuilder().denovoRef(0.0).denovoNonRef(0.0).create();

    //Check that adding new family members with 0 coverage doesn't change scores.
    final String[] base = {"", "", "AAAAAAACCCCCC", "AC"};
    final AbstractFamilyPosterior first = makeFamily(base); //, "", "", "C", "", "", "A", "AAACCC");
    final double basePosterior = first.getNonIdentityPosterior();
    for (int i = 1; i < 10; ++i) {
      final String[] next = new String[base.length + i];
      System.arraycopy(base, 0, next, 0, base.length);
      for (int pos = base.length; pos < next.length; ++pos) {
        next[pos] = "";
      }
      final AbstractFamilyPosterior incremental = makeFamily(next); //, "", "", "C", "", "", "A", "AAACCC");
      assertEquals(basePosterior, incremental.getNonIdentityPosterior(), 0.0000001);
    }
  }
  public void testFamilyScalingDenovo() throws Exception {
    mPriors = new GenomePriorParamsBuilder().create();

    // test case inspired by calls where all children are "de novo"
    // it isn't clear what the answer should be here.

    //Check that adding new family members with 0 coverage doesn't change scores.
    final String[] base = {"", "", "AAAAAAACCCCCC", "AC"};
    final AbstractFamilyPosterior first = makeFamily(base); //, "", "", "C", "", "", "A", "AAACCC");
    final double basePosterior = first.getNonIdentityPosterior();
    for (int i = 1; i < 10; ++i) {
      final String[] next = new String[base.length + i];
      System.arraycopy(base, 0, next, 0, base.length);
      for (int pos = base.length; pos < next.length; ++pos) {
        next[pos] = "";
      }
      final AbstractFamilyPosterior incremental = makeFamily(next); //, "", "", "C", "", "", "A", "AAACCC");
      assertEquals(basePosterior, incremental.getNonIdentityPosterior(), 0.0000001);
    }
  }
}
