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

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.reference.Ploidy;
import com.rtg.relation.Family;
import com.rtg.variant.GenomePriorParams;
import com.rtg.variant.bayes.Hypotheses;
import com.rtg.variant.bayes.ModelInterface;
import com.rtg.variant.bayes.multisample.HaploidDiploidHypotheses;
import com.rtg.variant.bayes.multisample.HypothesisScore;
import com.rtg.variant.bayes.snp.HypothesesNone;
import com.rtg.variant.bayes.snp.HypothesesPrior;

/**
 */
@TestClass(value = {"com.rtg.variant.bayes.multisample.family.FamilyPosteriorTest"})
public abstract class AbstractFamilyPosterior {

  protected final GenomePriorParams mParams;
  protected final ModelInterface<?> mFather;
  protected final ModelInterface<?> mMother;
  protected final List<ModelInterface<?>> mChildren;
  protected final List<HypothesisScore> mBest = new ArrayList<>();
  protected final int mHypothesesFatherSize;
  protected final int mHypothesesMotherSize;
  protected final int mReferenceHypothesis;
  protected final HypothesesPrior<?> mFatherHypotheses;
  protected final HypothesesPrior<?> mMotherHypotheses;
  protected final Ploidy mFatherPloidy;
  protected final Ploidy mMotherPloidy;
  protected final HaploidDiploidHypotheses<? extends HypothesesPrior<?>> mHypotheses;
  protected final double mLogDenovoRefPrior;
  protected final double mLogDenovoNonrefPrior;

  protected AbstractFamilyPosterior(Family family, GenomePriorParams params, List<ModelInterface<?>> models, HaploidDiploidHypotheses<? extends HypothesesPrior<?>> hypotheses) {
    final int[] ids = family.getSampleIds(); // Use IDs to pull selected family member individual models out of the model list
    mFather = models.get(ids[Family.FATHER_INDEX]);
    mMother = models.get(ids[Family.MOTHER_INDEX]);
    mHypotheses = hypotheses;
    mFatherHypotheses = hypotheses.get(mFather);
    mMotherHypotheses = hypotheses.get(mMother);
    if (mMother.reference() != mMotherHypotheses.reference()) {
      throw new IllegalArgumentException("All broken");
    }
    mHypothesesFatherSize = mFatherHypotheses.size();
    mHypothesesMotherSize = mMotherHypotheses.size();
    mReferenceHypothesis = mFatherHypotheses.reference(); // should be same for all family members
    assert mMotherHypotheses instanceof HypothesesNone || mReferenceHypothesis == mMotherHypotheses.reference();
    mFatherPloidy = mFatherHypotheses.ploidy();
    mMotherPloidy = mMotherHypotheses.ploidy();
    //assert mMother.hypotheses() == hypotheses; // not true for sex
    mChildren = new ArrayList<>(ids.length - Family.FIRST_CHILD_INDEX);
    for (int i = Family.FIRST_CHILD_INDEX; i < ids.length; ++i) {
      final ModelInterface<?> model = models.get(ids[i]);
      mChildren.add(model);
      final Hypotheses<?> hyp = model.hypotheses();
      assert hyp instanceof HypothesesNone || mReferenceHypothesis == hyp.reference();
    }
    mParams = params;
    mLogDenovoRefPrior = params.logDenovoRef();
    mLogDenovoNonrefPrior = params.logDenovoNonRef();
  }

  /**
   * @return true iff the family as a whole is not equal to reference.
   */
  public abstract boolean isInteresting();

  /**
   * @return ln of likelihood that this is a non-identity call.
   */
  public abstract double getNonIdentityPosterior();

  /**
   * @return the father's best category and posterior score
   */
  public HypothesisScore bestFather() {
    return mBest.get(Family.FATHER_INDEX);
  }

  /**
   * @return the mother's best category and posterior score
   */
  public HypothesisScore bestMother() {
    return mBest.get(Family.MOTHER_INDEX);
  }

  /**
   * @param i the index of the child
   * @return the child's best category and posterior score
   */
  public HypothesisScore bestChild(int i) {
    return mBest.get(i + Family.FIRST_CHILD_INDEX);
  }

  /**
   * @return the number of children (&ge; 1).
   */
  public int numberChildren() {
    return mChildren.size();
  }

}
