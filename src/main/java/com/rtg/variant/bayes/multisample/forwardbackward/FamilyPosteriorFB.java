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

package com.rtg.variant.bayes.multisample.forwardbackward;

import java.util.ArrayList;
import java.util.List;

import com.rtg.reference.Ploidy;
import com.rtg.relation.Family;
import com.rtg.util.MathUtils;
import com.rtg.variant.GenomePriorParams;
import com.rtg.variant.VariantSample;
import com.rtg.variant.bayes.ArrayGenotypeMeasure;
import com.rtg.variant.bayes.Code;
import com.rtg.variant.bayes.Factor;
import com.rtg.variant.bayes.ModelInterface;
import com.rtg.variant.bayes.multisample.HaploidDiploidHypotheses;
import com.rtg.variant.bayes.multisample.HypothesisScore;
import com.rtg.variant.bayes.multisample.family.AbstractFamilyPosterior;
import com.rtg.variant.bayes.multisample.family.MendelianAlleleProbability;
import com.rtg.variant.bayes.multisample.family.MendelianAlleleProbabilityFactory;
import com.rtg.variant.bayes.snp.HypothesesNone;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;

/**
 */
public class FamilyPosteriorFB extends AbstractFamilyPosterior {

  final PossibilityArithmetic mArithmetic;
  Family mFamily;
  final int[] mSampleIds;

  boolean mEqual;
  protected double mNonIdentity;
  protected final Code mMaximalCode;
  protected final MendelianAlleleProbabilityFactory mDenovoSubstitutionFactory;
  private Factor<?> mFatherMarginal;
  private Factor<?> mMotherMarginal;
  private Factor<?>[] mChildMarginal;

  FamilyPosteriorFB(Family family, GenomePriorParams params, List<ModelInterface<?>> models, HaploidDiploidHypotheses<?> hypotheses, MendelianAlleleProbabilityFactory denovoSubstitutionFactory) {
    super(family, params, models, hypotheses);
    mFamily = family;
    mSampleIds = mFamily.getSampleIds();
    mArithmetic = mFatherHypotheses.arithmetic();
    assert mArithmetic == mMotherHypotheses.arithmetic() : mArithmetic.getClass().getName() + " != " + mMotherHypotheses.arithmetic().getClass().getName();
    final double zero = mArithmetic.zero();
    mNonIdentity = zero;


    mMaximalCode = mFatherHypotheses.code().size() > mMotherHypotheses.code().size() ? mFatherHypotheses.code() : mMotherHypotheses.code();
 //   assert !hypotheses.isDefault(); // we're ignoring this for now?!?
    mDenovoSubstitutionFactory = denovoSubstitutionFactory;
  }

  /**
   * Testing context only
   * @return initial values for A when no A calculations have been done on samples
   */
  Factor<?>[] initialA() {
    final Factor<?>[] as = new Factor<?>[mFamily.size()];
    as[Family.FATHER_INDEX] = CommonFormulas.createMutableFactor(mFatherHypotheses);
    as[Family.MOTHER_INDEX] = CommonFormulas.createMutableFactor(mMotherHypotheses);
    for (int i = 0; i < mChildren.size(); ++i) {
      as[Family.FIRST_CHILD_INDEX + i] = null; //CommonFormulas.hypothesesToVector(mChildren.get(i).hypotheses()); // possibly not needed
    }
    return as;
  }

  /**
   * Testing context only, assumes single generation
   * @return initial values for B when no B calculations have been done on samples
   */
  BContainer[] initialB() {
    final ArrayList<ModelInterface<?>> models = new ArrayList<>(mFamily.size());
    final int[] bsize = new int[mFamily.size()];
    bsize[Family.FATHER_INDEX] = mFamily.getFatherDistinctMates();
    bsize[Family.MOTHER_INDEX] = mFamily.getMotherDistinctMates();
    models.add(mFather);
    models.add(mMother);
    for (int i = 0; i < mChildren.size(); ++i) {
      models.add(mChildren.get(i));
      bsize[Family.FIRST_CHILD_INDEX + i] = 1;
    }
    return CommonFormulas.initialB(models, bsize);
  }

  Factor<?>[] computeChildAs(Factor<?>[] as, BContainer[] bs) {
    return new ACalculator(as, bs).calculateA();
  }

  Factor<?>[] computeParentBs(Factor<?>[] as, BContainer[] bs) {
    return new BCalculator(as, bs).calculateB();
  }

  void computeMarginals(Factor<?>[] as, BContainer[] bs) {
    final Factor<?> fatherModelVector = CommonFormulas.createMutableFactor(mFather);
    mFatherMarginal = CommonFormulas.computeQ(as[mSampleIds[Family.FATHER_INDEX]], fatherModelVector, bs[mSampleIds[Family.FATHER_INDEX]]);
    mMotherMarginal = CommonFormulas.computeQ(as[mSampleIds[Family.MOTHER_INDEX]], CommonFormulas.createMutableFactor(mMother), bs[mSampleIds[Family.MOTHER_INDEX]]);
    mChildMarginal = new Factor<?>[mChildren.size()];
    for (int i = 0; i < mChildMarginal.length; ++i) {
      mChildMarginal[i] = CommonFormulas.computeQ(as[mSampleIds[Family.FIRST_CHILD_INDEX + i]], CommonFormulas.createMutableFactor(mChildren.get(i)), bs[mSampleIds[Family.FIRST_CHILD_INDEX + i]]);
    }
  }

  void addMarginals(FamilyPosteriorFB other) {
    mFatherMarginal = addVector(mFatherMarginal, other.mFatherMarginal);
    mMotherMarginal = addVector(mMotherMarginal, other.mMotherMarginal);

    if (mChildMarginal == null) {
      mChildMarginal = new Factor<?>[mChildren.size()];
    }
    for (int i = 0; i < mChildren.size(); ++i) {
      mChildMarginal[i] = addVector(mChildMarginal[i], other.mChildMarginal[i]);
    }
  }

  Factor<?> addVector(Factor<?> a, Factor<?> b) {
    if (a == null) {
      return b;
    }
    final PossibilityArithmetic arithmetic = a.arithmetic();
    final MutableFactor<?> mutableFactor = new MutableFactor<>(a.hypotheses(), arithmetic, a.size());
    for (int i = 0; i < a.size(); ++i) {
      mutableFactor.set(i, arithmetic.add(a.p(i), arithmetic.poss2Poss(b.p(i), b.arithmetic())));
    }
    return mutableFactor;
  }
  double[] factorToDoubles(Factor<?> fact) {
    final double[] res = new double[fact.size()];
    for  (int i = 0; i < fact.size(); ++i) {
      res[i] = fact.p(i);
    }
    return res;
  }

  /**
   * Create the best hypothesis scores for the family
   * @param denovoPosteriors de novo scores for the children
   */
  public void findBest(double[] denovoPosteriors) {
    final HypothesisScore father = mFatherMarginal.hypotheses().size() == 0 ? null : new HypothesisScore(new ArrayGenotypeMeasure(mArithmetic, factorToDoubles(mFatherMarginal), mFatherHypotheses));
    mBest.add(father);
    final HypothesisScore mother = mMotherMarginal.hypotheses().size() == 0 ? null : new HypothesisScore(new ArrayGenotypeMeasure(mArithmetic, factorToDoubles(mMotherMarginal), mMotherHypotheses));
    mBest.add(mother);
    boolean anyDenovo = false;
    for (int i = 0; i < mChildren.size(); ++i) {
      final HypothesisScore child = mChildMarginal[i].hypotheses().size() == 0 ? null : new HypothesisScore(new ArrayGenotypeMeasure(mArithmetic, factorToDoubles(mChildMarginal[i]), mChildren.get(i).hypotheses()));
      mBest.add(child);
      if (child != null) {
        final Ploidy childPloidy = mHypotheses.get(mChildren.get(i)).ploidy();
        final int fatherHyp = father == null ? -1 : father.hypothesis();
        final int motherHyp = mother == null ? -1 : mother.hypothesis();
        final int reference = mHypotheses.get(mChildren.get(i)).reference();
        final MendelianAlleleProbability al = MendelianAlleleProbabilityFactory.COMBINED.getMendelianAlleleProbability(mFatherPloidy, mMotherPloidy, childPloidy, mLogDenovoRefPrior, mLogDenovoNonrefPrior, reference);
        if (al.isDenovo(mMaximalCode, fatherHyp, motherHyp, child.hypothesis())) {
          child.setDenovo(VariantSample.DeNovoStatus.IS_DE_NOVO);
          anyDenovo = true;
        }
        child.setDeNovoPosterior(denovoPosteriors[i]);
      }
    }
    double sumNips = 0.0;
    boolean equal = true;
    for (int i = 0; i < mBest.size(); ++i) {
      final HypothesisScore best = mBest.get(i);
      if (best != null) {
        equal &= mReferenceHypothesis == best.hypothesis();
        final double x = best.nonIdentityPosterior();
        sumNips += MathUtils.logExpPlus1(x); // Math.log(Math.exp(nip) + 1), nip = -ip for ref= callsa
        if (anyDenovo) {
          if (i >= Family.FIRST_CHILD_INDEX && best.isDeNovo() == VariantSample.DeNovoStatus.UNSPECIFIED) { //if any child is denovo go through the family and set denovo flag for children
            best.setDenovo(VariantSample.DeNovoStatus.NOT_DE_NOVO);
          }
        }
      }
    }


    mNonIdentity = MathUtils.logExpMinus1(sumNips);
    mEqual = equal;
  }

  /**
   * Testing context only
   * Compute family marginals assuming priors have already corrected for number of alleles. E.g. via site-specific priors or via EM optimization
   */
  void computeAlleles() {
    final Factor<?>[] as = initialA();
    final BContainer[] bs = initialB();
    computeAlleles(as, bs);
  }

  private void computeAlleles(Factor<?>[] as, BContainer[] bs) {
    final Factor<?>[] childAs = computeChildAs(as, bs);
    System.arraycopy(childAs, 0, as, 2, childAs.length);
    final Factor<?>[] parentBs = computeParentBs(as, bs);
    bs[Family.FATHER_INDEX].setB(mFamily.getFatherFamilyId(), parentBs[Family.FATHER_INDEX]);
    bs[Family.MOTHER_INDEX].setB(mFamily.getMotherFamilyId(), parentBs[Family.MOTHER_INDEX]);
    computeMarginals(as, bs);
    findBest(new double[mChildren.size()]);
  }

  /**
   * Return the sum of the specified child's marginals
   * @param child a child index
   * @return the sum of all marginals for the child
   */
  public double childMarginalSum(int child) {
    final Factor<?> vector = mChildMarginal[child];
    return marginalSum(vector);
  }

  private double marginalSum(Factor<?> hypothesesVector) {
    final PossibilityArithmetic arithmetic = hypothesesVector.arithmetic();
    double sum = arithmetic.zero();
    for (int i = 0; i < hypothesesVector.size(); ++i) {
      sum = arithmetic.add(sum, hypothesesVector.p(i));
    }
    return arithmetic.poss2Ln(sum);
  }

  @Override
  public boolean isInteresting() { // TODO is this correct surely low posterior reference calls are also interesting?
    return !mEqual;
  }

  @Override
  public double getNonIdentityPosterior() {
    return mNonIdentity;
  }

  /**
   * @return the father's best category and posterior score
   */
  @Override
  public HypothesisScore bestFather() {
    return mBest.get(Family.FATHER_INDEX);
  }

  /**
   * @return the mother's best category and posterior score
   */
  @Override
  public HypothesisScore bestMother() {
    return mBest.get(Family.MOTHER_INDEX);
  }
  /**
   * @param i the index of the child
   * @return the child's best category and posterior score
   */
  @Override
  public HypothesisScore bestChild(int i) {
    return mBest.get(i + Family.FIRST_CHILD_INDEX);
  }

  private class InnerCalculator {
    final Factor<?> mFatherA;
    final Factor<?> mFatherS;
    final Factor<?> mFatherE;
    final Factor<?> mMotherA;
    final Factor<?> mMotherS;
    final Factor<?> mMotherE;
    final Factor<?>[] mD;
    final double[][][] mC;

    InnerCalculator(Factor<?>[] as, BContainer[] bs) {
      mFatherA = as[mSampleIds[Family.FATHER_INDEX]];
      mFatherS = CommonFormulas.createMutableFactor(mFather);
      mFatherE = CommonFormulas.forwardEGeneral(mFamily.getFatherFamilyId(), mFatherA, mFatherS, bs[mSampleIds[Family.FATHER_INDEX]]);

      mMotherA = as[mSampleIds[Family.MOTHER_INDEX]];
      mMotherS = CommonFormulas.createMutableFactor(mMother);
      mMotherE = CommonFormulas.forwardEGeneral(mFamily.getMotherFamilyId(), mMotherA, mMotherS, bs[mSampleIds[Family.MOTHER_INDEX]]);

      final int childSize = mChildren.size();
      mD = new Factor<?>[childSize];
      mC = new double[childSize][][];
      for (int a = 0; a < childSize; ++a) {
        final ModelInterface<?> model = mChildren.get(a);
        final Factor<?> s = CommonFormulas.createMutableFactor(model);
        final BContainer b =  bs[mSampleIds[Family.FIRST_CHILD_INDEX + a]];
        mD[a] = CommonFormulas.computeD(s, b);
        // mDenovoNonrefPrior, mLogDenovoNonrefPrior
        final MendelianAlleleProbability m = MendelianAlleleProbabilityFactory.COMBINED.getMendelianAlleleProbability(mFatherPloidy, mMotherPloidy, model.hypotheses().ploidy(), mLogDenovoRefPrior, mLogDenovoNonrefPrior, model.reference());
        mC[a] = CommonFormulas.backwardC(mFatherA.size(), mMotherA.size(), CommonFormulas.maxCode(mFatherHypotheses, mMotherHypotheses) , mD[a], m);
      }
    }
  }

  private class ACalculator extends InnerCalculator {
    ACalculator(Factor<?>[] as, BContainer[] bs) {
      super(as, bs);
    }

    Factor<?>[] calculateA() {
      final int childSize = mChildren.size();
      final Factor<?>[] ret = new Factor<?>[childSize];
      for (int a = 0; a < childSize; ++a) {
        final ModelInterface<?> model = mChildren.get(a);
        final MendelianAlleleProbability m = mDenovoSubstitutionFactory.getMendelianAlleleProbability(mFatherPloidy, mMotherPloidy, model.hypotheses().ploidy(), mLogDenovoRefPrior, mLogDenovoNonrefPrior, model.reference());
        final Factor<?> av = CommonFormulas.forwardA(mFatherE, mMotherE, model.hypotheses(), a, mC, m);
        ret[a] = av;
        //mChildMarginal[a] = CommonFormulas.dot(av, d[a]);
      }
      return ret;
    }
  }
  private class BCalculator extends InnerCalculator {
    final BContainer[] mBs;
    BCalculator(Factor<?>[] as, BContainer[] bs) {
      super(as, bs);
      mBs = bs;
    }

    Factor<?>[] calculateB() {
      final Factor<?> fatherB;
      if (!(mFatherHypotheses instanceof HypothesesNone)) {
        fatherB = CommonFormulas.backwardB(mFatherE, mMotherE, mC, true);
      } else {
        fatherB = mBs[mSampleIds[Family.FATHER_INDEX]].getB(mFamily.getFatherFamilyId());
      }
      final Factor<?> motherB;
      if (!(mMotherHypotheses instanceof  HypothesesNone)) {
        motherB = CommonFormulas.backwardB(mMotherE, mFatherE, mC, false);
      } else {
        motherB = mBs[mSampleIds[Family.MOTHER_INDEX]].getB(mFamily.getMotherFamilyId());
      }
      return new Factor<?>[]{fatherB, motherB};
    }
  }
}
