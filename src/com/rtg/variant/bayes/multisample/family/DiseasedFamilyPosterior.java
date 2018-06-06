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
import java.util.Arrays;
import java.util.List;
import java.util.Stack;

import com.rtg.relation.Family;
import com.rtg.variant.GenomePriorParams;
import com.rtg.variant.bayes.ArrayGenotypeMeasure;
import com.rtg.variant.bayes.Code;
import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.Hypotheses;
import com.rtg.variant.bayes.ModelInterface;
import com.rtg.variant.bayes.NoNonIdentityMeasure;
import com.rtg.variant.bayes.multisample.HypothesisScore;
import com.rtg.variant.bayes.multisample.Utils;
import com.rtg.variant.util.VariantUtils;
import com.rtg.variant.util.arithmetic.LogApproximatePossibility;

/**
 * Calculates the posterior probability of an allele accounting for a disease.
 * We make some assumptions:
 * <ul>
 * <li> the allele is dominant</li>
 * <li> that only one parent has the disease </li>
 * <li> the disease is rare </li>
 * <li> the reference sequence doesn't have the disease </li>
 * <li> the disease status of each family member is known </li>
 * </ul>
 * XXX We aren't actually computing a proper GQ for this caller. We're computing P(H_a|H_d, E) instead of P(H_a|E)
 * There is no Non Identity posterior score either
 */
public class DiseasedFamilyPosterior {

  private HypothesisScore marginal(double[][] marginals, int disease, Hypotheses<?> hypotheses) {
    final double[] result = new double[marginals.length];
    for (int i = 0; i < marginals.length; ++i) {
      result[i] = marginals[i][disease];
    }
    return new HypothesisScore(new NoNonIdentityMeasure(new ArrayGenotypeMeasure(LogApproximatePossibility.SINGLETON, result, hypotheses)));
  }

  final GenomePriorParams mParams;
  final ModelInterface<?> mFather;
  final ModelInterface<?> mMother;
  final List<ModelInterface<?>> mChildren;
  final double[][] mFatherMarginal;
  final double[][] mMotherMarginal;
  final double[] mDiseaseMarginal;
  final List<double[][]> mChildMarginal;
  final Family mFamily;

  private final List<HypothesisScore> mBest = new ArrayList<>();
  private HypothesisScore mBestDisease = null;
  protected final Hypotheses<?> mHypotheses;
  private final int mHypothesesSize;
  protected final Code mCode;
  protected final HypothesesDisease mDiseaseHypotheses;

  //TODO add flag for recessive vs dominant
  DiseasedFamilyPosterior(GenomePriorParams params, Family family, Hypotheses<?> commonHypotheses, HypothesesDisease diseaseHypotheses, List<ModelInterface<?>> models) {
    mHypotheses = commonHypotheses;
    mCode = mHypotheses.code();
    mHypothesesSize = mHypotheses.size();
    final int[] ids = family.getSampleIds();
    mFather = models.get(ids[Family.FATHER_INDEX]);
    mMother = models.get(ids[Family.MOTHER_INDEX]);
    assert mMother.hypotheses() == mHypotheses;
    assert mFather.hypotheses() == mHypotheses;
    mChildren = new ArrayList<>(ids.length - Family.FIRST_CHILD_INDEX);
    for (int i = Family.FIRST_CHILD_INDEX; i < ids.length; ++i) {
      final ModelInterface<?> model = models.get(ids[i]);
      mChildren.add(model);
      assert model.hypotheses() == mHypotheses;
    }
    mParams = params;
    mFamily = family;

    mDiseaseHypotheses = diseaseHypotheses;
    final Description diseaseDescription = mDiseaseHypotheses.description();
    final int diseaseSize = diseaseDescription.size();
    mFatherMarginal = new double[mHypothesesSize][diseaseSize];
    mMotherMarginal = new double[mHypothesesSize][diseaseSize];
    mDiseaseMarginal = new double[diseaseSize];
    mChildMarginal = new ArrayList<>();
    for (final double[] x : mFatherMarginal) {
      Arrays.fill(x, Double.NEGATIVE_INFINITY);
    }
    for (final double[] x : mMotherMarginal) {
      Arrays.fill(x, Double.NEGATIVE_INFINITY);
    }
    Arrays.fill(mDiseaseMarginal, Double.NEGATIVE_INFINITY);
    for (int k = 0; k < mChildren.size(); ++k) {
      final double[][] childMarginals = new double[mHypothesesSize][diseaseSize];
      for (final double[] x : childMarginals) {
        Arrays.fill(x, Double.NEGATIVE_INFINITY);
      }
      mChildMarginal.add(childMarginals);
    }

  }

  void process() {
    for (int i = 0; i < mHypothesesSize; ++i) {
      for (int j = 0; j < mHypothesesSize; ++j) {
        double marginal = 0.0;
        // p dagger
        final double generatedPrior =  AlleleSetProbabilityDiploid.SINGLETON.probabilityLn(mCode, mHypotheses.reference(), i, j);
        marginal += mFather.posteriorLn0(i);
        marginal += mMother.posteriorLn0(j);
        marginal += generatedPrior;
        calculateMarginals(marginal, i, j);
      }
    }
    mBestDisease = new HypothesisScore(new ArrayGenotypeMeasure(LogApproximatePossibility.SINGLETON, mDiseaseMarginal, mDiseaseHypotheses));

    // Marginals given we have chosen the disease explanation
    final int diseaseIndex = mBestDisease.hypothesis();
    mBest.add(marginal(mFatherMarginal, diseaseIndex, mFather.hypotheses()));
    mBest.add(marginal(mMotherMarginal, diseaseIndex, mMother.hypotheses()));
    for (int i = 0; i < mChildMarginal.size(); ++i) {
      mBest.add(marginal(mChildMarginal.get(i), diseaseIndex, mChildren.get(i).hypotheses()));
    }
  }

  /**
   * @param marginal accumulated score so far
   * @param father current hypothesis for the father
   * @param mother current hypothesis for the mother
   */
  void calculateMarginals(double marginal, int father, int mother) {
    calculateMarginalsSlowly(marginal, father, mother, new Stack<>(), 0);
  }

  private void calculateMarginalsSlowly(double marginal, int father, int mother, Stack<Integer> childCats, int childCount) {
    if (childCount == mChildren.size()) {
      updateMarginals(marginal, father, mother, childCats);
    } else {
      final ModelInterface<?> child = mChildren.get(childCount);
      for (int i = 0; i < mHypothesesSize; ++i) {
        final double prior = MendelianAlleleProbabilityDiploid.SINGLETON.probabilityLn(mCode, father, mother, i);
        if (Double.isInfinite(prior)) {
          continue;
        }
        final double childMarginal = marginal + child.posteriorLn0(i) + prior;
        childCats.push(i);
        calculateMarginalsSlowly(childMarginal, father, mother, childCats, childCount + 1);
        childCats.pop();
      }
    }
  }

  /**
   * @param b accumulated score so far
   * @param father current hypothesis for the father
   * @param mother current hypothesis for the mother
   * @param childHyps hypotheses for the children
   */
  private void updateMarginals(final double b, int father, int mother, Stack<Integer> childHyps) {
    final int[] currentHyp = new int[childHyps.size() + Family.FIRST_CHILD_INDEX];
    currentHyp[Family.FATHER_INDEX] = father;
    currentHyp[Family.MOTHER_INDEX] = mother;
    for (int i = 0; i < childHyps.size(); ++i) {
      currentHyp[i + Family.FIRST_CHILD_INDEX] = childHyps.get(i);
    }
    final Disease disease = new Disease(mFamily, mHypotheses, currentHyp);

    // Either there is 1 explanation or NO_DISEASE at this position
    final int explanation = disease.explanation(mHypotheses.reference());
    final double diseasePrior = mDiseaseHypotheses.prior(explanation);

    // p prime
    final int alleleCountExcl = Utils.numberAllelesExclude(mHypotheses.description().size(), mCode, explanation - 1, father, mother, mHypotheses.reference());
    //System.err.println("father=" + mHypotheses.name(father) + " mother=" + mHypotheses.name(mother) + " expl=" + mDiseaseHypotheses.name(explanation) + " cnt= " + alleleCountExcl);
    final double pPrime = mParams.getAlleleFrequencyLnProbability(alleleCountExcl);
    final double finalMarginal = b + pPrime + diseasePrior;
    //update marginals.
    mFatherMarginal[father][explanation] = VariantUtils.logSumApproximation(mFatherMarginal[father][explanation], finalMarginal);
    mMotherMarginal[mother][explanation] = VariantUtils.logSumApproximation(mMotherMarginal[mother][explanation], finalMarginal);
    //    System.err.println("Xfather=" + father + " mother=" + mother + " fm=" + finalMarginal);
    mDiseaseMarginal[explanation] = VariantUtils.logSumApproximation(mDiseaseMarginal[explanation], finalMarginal);
    //    System.err.println("d=" + d + " father=" + father + " mother=" + mother + " b=" + b + " pPrime=" + pPrime + " q=" + q + " set=" + alleleSet + " dm[d]=" + mDiseaseMarginal[d] + " nt=" + refNt);
    for (int i = 0; i < childHyps.size(); ++i) {
      final double[][] marginals = mChildMarginal.get(i);
      final int cat = childHyps.get(i);
      marginals[cat][explanation] = VariantUtils.logSumApproximation(marginals[cat][explanation], finalMarginal);
    }
  }

  public boolean isInteresting() {
    return mBestDisease.hypothesis() != 0;
  }

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
   * @return the best explanation of the disease and posterior score
   */
  public HypothesisScore bestDisease() {
    return mBestDisease;
  }

  /**
   * @param i the index of the child
   * @return the child's best category and posterior score
   */
  public HypothesisScore bestChild(int i) {
    return mBest.get(Family.FIRST_CHILD_INDEX + i);
  }

  double anyDiseasePosteriorRatio() {
    double explanation = Double.NEGATIVE_INFINITY;
    // position 0 is the no-disease case
    for (int k = 1; k < mDiseaseMarginal.length; ++k) {
      explanation = VariantUtils.logSumApproximation(explanation, mDiseaseMarginal[k]);
    }
    return explanation - mDiseaseMarginal[0];
  }
}
