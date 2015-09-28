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

import static com.rtg.variant.VariantSample.DeNovoStatus;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Stack;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.reference.Ploidy;
import com.rtg.relation.Family;
import com.rtg.util.QuickSort;
import com.rtg.util.QuickSortDoubleIntProxy;
import com.rtg.util.diagnostic.SpyHistogram;
import com.rtg.variant.GenomePriorParams;
import com.rtg.variant.bayes.AlleleStatistics;
import com.rtg.variant.bayes.ArrayGenotypeMeasure;
import com.rtg.variant.bayes.Code;
import com.rtg.variant.bayes.Hypotheses;
import com.rtg.variant.bayes.ModelInterface;
import com.rtg.variant.bayes.multisample.HaploidDiploidHypotheses;
import com.rtg.variant.bayes.multisample.HypothesisScore;
import com.rtg.variant.bayes.multisample.Utils;
import com.rtg.variant.bayes.snp.HypothesesPrior;
import com.rtg.variant.util.VariantUtils;
import com.rtg.variant.util.arithmetic.LogApproximatePossibility;

/**
 * Computation of posteriors for a family.
 */
@TestClass(value = {"com.rtg.variant.bayes.multisample.family.FamilyPosteriorTest", "com.rtg.variant.bayes.multisample.family.FastFamilyPosteriorTest"})
public class FamilyPosterior extends AbstractFamilyPosterior {

  protected static final int MISSING_HYPOTHESIS = -1;

  private static final double TERM_THRESHOLD = 10.0;
  static final boolean TERMINATE = true;
  static final boolean SQUARE_FIRST = true;
  static final boolean ENABLE_SORTED_HYPOTHESES = true;
  private static final boolean FIRST_HYP_FULLY = true;

  double[] mFatherMarginal;
  double[] mMotherMarginal;
  List<double[]> mChildMarginal;
  double[] mChildDenovoMarginal = null;
  double[] mChildNonDenovoMarginal = null;
  protected final Code mMaximalCode;
  protected final boolean mEqual;
  protected double mNonIdentity = Double.NEGATIVE_INFINITY;
  protected double mIdentity = Double.NEGATIVE_INFINITY;
  private final double mContraryProbabilityLn;

  FamilyPosterior(Family family, GenomePriorParams params, List<ModelInterface<?>> models, HaploidDiploidHypotheses<?> hypotheses) {
    super(family, params, models, hypotheses);
    mMaximalCode = mFatherHypotheses.code().size() > mMotherHypotheses.code().size() ? mFatherHypotheses.code() : mMotherHypotheses.code();
    mContraryProbabilityLn = Math.log(params.contraryProbability());

    if (mHypothesesFatherSize == 0 && mHypothesesMotherSize == 0) {
      mBest.add(null); // father
      mBest.add(null); // mother
      for (int i = 0; i < mChildren.size(); i++) {
        mBest.add(null);
      }
      mEqual = true;
    } else {
      boolean anyDenovo = false;
      // Because calculating a de novo marginal is slow (observed 40% slowdown)
      // we only want to do it in the case that we actually think there may be
      // a de novo mutation. These are rare so we can afford the complete
      // recompute
      List<HypothesisScore> bestList;
      do {
        bestList = new ArrayList<>();
        mFatherMarginal = createMarginalArray(mHypothesesFatherSize);
        mMotherMarginal = createMarginalArray(mHypothesesMotherSize);
        mChildMarginal = new ArrayList<>();
        if (anyDenovo) {
          mChildDenovoMarginal = createMarginalArray(mChildren.size());
          mChildNonDenovoMarginal = createMarginalArray(mChildren.size());
        }

        for (ModelInterface<?> child : mChildren) {
          final int size = child.size();
          mChildMarginal.add(createMarginalArray(size));
        }

        if (ENABLE_SORTED_HYPOTHESES && mFatherHypotheses.size() * mMotherHypotheses.size() >= 100) {
          computeSuperAlleles(hypotheses.isDefault());
        } else {
          computeAlleles(hypotheses.isDefault());
        }
        final HypothesisScore father = mFatherMarginal.length == 0 ? null : new HypothesisScore(new ArrayGenotypeMeasure(LogApproximatePossibility.SINGLETON, mFatherMarginal, mFather.hypotheses()));
        bestList.add(father);
        final HypothesisScore mother = mMotherMarginal.length == 0 ? null : new HypothesisScore(new ArrayGenotypeMeasure(LogApproximatePossibility.SINGLETON, mMotherMarginal, mMother.hypotheses()));
        bestList.add(mother);
        for (int i = 0; i < mChildren.size(); i++) {
          final ModelInterface<?> childModel = mChildren.get(i);
          final HypothesisScore best = mChildMarginal.get(i).length == 0 ? null
            : new HypothesisScore(new ArrayGenotypeMeasure(LogApproximatePossibility.SINGLETON, mChildMarginal.get(i), childModel.hypotheses()));
          if (mChildDenovoMarginal != null && best != null) {
            markDenovo(best, mChildDenovoMarginal[i], mChildNonDenovoMarginal[i]);
          }
          if (isDenovo(bestHypothesis(father), bestHypothesis(mother), i, bestHypothesis(best))) {
            anyDenovo = true;
          }
          bestList.add(best);
        }
      } while (anyDenovo && mChildDenovoMarginal == null);
      mBest.addAll(bestList);
      boolean equal = true;
      for (int i = 0; i < mBest.size(); i++) {
        final HypothesisScore best = mBest.get(i);
        if (best != null) {
          equal &= mReferenceHypothesis == best.hypothesis();
          if (anyDenovo) {
            if (i >= Family.FIRST_CHILD_INDEX && best.isDeNovo() == DeNovoStatus.UNSPECIFIED) { //if any child is denovo go through the family and set denovo flag for children
              best.setDenovo(DeNovoStatus.NOT_DE_NOVO);
            }
          }
        }
      }
      mEqual = equal;
    }
  }


  static int bestHypothesis(HypothesisScore best) {
    return best == null ? MISSING_HYPOTHESIS : best.hypothesis();
  }

  private static double[] createMarginalArray(int size) {
    final double[] childMarginals = new double[size];
    Arrays.fill(childMarginals, Double.NEGATIVE_INFINITY);
    return childMarginals;
  }

  @Override
  public boolean isInteresting() {
    return !mEqual;
  }

  /**
   * Compute family marginals
   * @param isDefault if true, calculate default priors, if false, assumes priors have already
   * corrected for number of alleles. E.g. via site-specific priors or via EM optimization
   */
  private void computeAlleles(boolean isDefault) {
    final int haploidSize;
    final int fatherIterations;
    final AlleleProbability ap = isDefault ? getAlleleProbability(mFatherPloidy, mMotherPloidy) : null;
    final boolean fatherNone = mFatherPloidy == Ploidy.NONE;
    if (fatherNone) {
      fatherIterations = 1;
      haploidSize = mMotherHypotheses.description().size();
    } else {
      fatherIterations = mHypothesesFatherSize;
      haploidSize = mFatherHypotheses.description().size();
    }
    final boolean motherNone = mMotherPloidy == Ploidy.NONE;
    final int motherIterations = motherNone ? 1 : mMotherHypotheses.size();
    for (int i = 0; i < fatherIterations; i++) {
      final double fatherpost;
      if (fatherNone) {
        fatherpost = 0;
      } else {
        fatherpost = mFather.posteriorLn0(i) + (isDefault ? 0 : mFatherHypotheses.arithmetic().poss2Ln(mFatherHypotheses.p(i)));
      }
      final int fatherHypothesis =  fatherNone ? MISSING_HYPOTHESIS : i;
      for (int j = 0; j < motherIterations; j++) {
        final double marginal;
        final double motherpost;
        final int motherHypothesis =  motherNone ? MISSING_HYPOTHESIS : j;
        if (motherNone) {
          motherpost = 0;
        } else {
          motherpost = mMother.posteriorLn0(j) + (isDefault ? 0 : mMotherHypotheses.arithmetic().poss2Ln(mMotherHypotheses.p(j)));
        }
        if (ap != null) {
          // p dagger in theory doc
          final double generatedPrior = ap.probabilityLn(mMaximalCode, mReferenceHypothesis, fatherHypothesis, motherHypothesis);
          // p prime in theory doc
          final int alleles = Utils.numberAlleles(haploidSize, mMaximalCode, mReferenceHypothesis, fatherHypothesis, motherHypothesis);
          marginal =
              -BinomialSpecial.logBinomial(haploidSize, alleles)
              + mParams.getAlleleFrequencyLnProbability(alleles)
              + fatherpost
              + motherpost
              + generatedPrior;
        } else {
          marginal = fatherpost + motherpost;
        }
//        System.err.println(fatherHypothesis + " " + motherHypothesis + " " + marginal);
        calculateMarginals(marginal, fatherHypothesis, motherHypothesis);
      }
    }
  }

  private QuickSortDoubleIntProxy makeSortedHypLookup(HypothesesPrior<?> hypotheses, ModelInterface<?> model, boolean isDefault) {
    final double[] scores = new double[hypotheses.size()];
    final int[] hypLookup = new int[hypotheses.size()];
    for (int i = 0; i < scores.length; i++) {
      scores[i] = model.posteriorLn0(i) + (isDefault ? 0 : hypotheses.arithmetic().poss2Ln(hypotheses.p(i)));
      hypLookup[i] = i;
    }
    final QuickSortDoubleIntProxy proxy = new QuickSortDoubleIntProxy(scores, hypLookup);
    QuickSort.sort(proxy);
    return proxy;
  }

  /**
   * Version of compute alleles that terminates early when it thinks it has enough precision.
   *
   * Compute marginals by considering hypotheses in order of most likely parent combinations then
   * estimating an upper bound on the un computed part of the calculation. If the upper bound is
   * sufficiently small relative to the computed part then we can stop work.
   *
   * Compute family marginals applying corrections to reduce the likelihood of multi-allelic sites
   * @param isDefault if true, calculate default priors, if false, assumes priors have already
   * corrected for number of alleles. E.g. via site-specific priors or via EM optimization
   */
  private void computeSuperAlleles(boolean isDefault) {
    final AlleleProbability ap = getAlleleProbability(mFatherPloidy, mMotherPloidy);
    assert mFatherPloidy != Ploidy.NONE || mMotherPloidy != Ploidy.NONE;
    final int haploidSize;
    final double[] fatherPost;
    final int[] fatherHyp;
    final double[] motherPost;
    final int[] motherHyp;
    if (mFatherPloidy == Ploidy.NONE) {
      haploidSize = mMotherHypotheses.description().size();
      fatherPost = new double[] {0};
      fatherHyp = new int[] {MISSING_HYPOTHESIS};
    } else {
      haploidSize = mFatherHypotheses.description().size();
      final QuickSortDoubleIntProxy proxy = makeSortedHypLookup(mFatherHypotheses, mFather, isDefault);
      fatherPost = proxy.getVals();
      fatherHyp = proxy.getPairArray();
    }
    if (mMotherPloidy == Ploidy.NONE) {
      motherPost = new double[] {0};
      motherHyp = new int[] {MISSING_HYPOTHESIS};
    } else {
      final QuickSortDoubleIntProxy proxy = makeSortedHypLookup(mMotherHypotheses, mMother, isDefault);
      motherPost = proxy.getVals();
      motherHyp = proxy.getPairArray();
    }

    final double[] fatherCum = makeCumulative(fatherPost);
    final double[] motherCum = makeCumulative(motherPost);

    final int[] nextFather = new int[fatherPost.length];
    final double[] faHypRem = Arrays.copyOf(fatherPost, fatherPost.length);
    final double[] moHypRem = Arrays.copyOf(motherPost, motherPost.length);

    int motherReferenceIndex = 0;
    for (int i = 0; i < motherHyp.length; i++) {
      if (motherHyp[i] == mMotherHypotheses.reference()) {
        motherReferenceIndex = i;
      }
    }
    int fatherReferenceIndex = 0;
    for (int i = 0; i < fatherHyp.length; i++) {
      if (fatherHyp[i] == mFatherHypotheses.reference()) {
        fatherReferenceIndex = i;
      }
    }
    int iterations = 0;

    if (FIRST_HYP_FULLY) {
      for (int i = 0; i < fatherHyp.length; i++) {
        final int j = 0;
        computeParentMarginal(isDefault, ap, haploidSize, fatherPost[i], motherPost[j], fatherHyp[i], motherHyp[j]);
        faHypRem[i] = j == motherCum.length - 1 ? Double.NEGATIVE_INFINITY : fatherPost[i] + motherCum[j + 1];
        moHypRem[j] = i == fatherCum.length - 1 ? Double.NEGATIVE_INFINITY : motherPost[j] + fatherCum[i + 1];
        nextFather[i]++;
        iterations++;
      }
      for (int j = 1; j < motherHyp.length; j++) {
        final int i = 0;
        computeParentMarginal(isDefault, ap, haploidSize, fatherPost[i], motherPost[j], fatherHyp[i], motherHyp[j]);
        faHypRem[i] = j == motherCum.length - 1 ? Double.NEGATIVE_INFINITY : fatherPost[i] + motherCum[j + 1];
        moHypRem[j] = i == fatherCum.length - 1 ? Double.NEGATIVE_INFINITY : motherPost[j] + fatherCum[i + 1];
        nextFather[i]++;
        iterations++;
      }
    }
    if (SQUARE_FIRST) {
      // Optimization.
      // Compute a square in the most likely portion of the grid first without doing early termination checks
      final int fatherLimit = Math.min(Math.max(fatherReferenceIndex, 1), fatherPost.length); //XXX
      final int motherLimit = Math.min(Math.max(motherReferenceIndex, 1), motherPost.length); //XXX
      for (int i = FIRST_HYP_FULLY ? 1 : 0; i <= fatherLimit; i++) {
        for (int j = FIRST_HYP_FULLY ? 1 : 0; j <= motherLimit; j++) {
          computeParentMarginal(isDefault, ap, haploidSize, fatherPost[i], motherPost[j], fatherHyp[i], motherHyp[j]);
          faHypRem[i] = j == motherCum.length - 1 ? Double.NEGATIVE_INFINITY : fatherPost[i] + motherCum[j + 1];
          moHypRem[j] = i == fatherCum.length - 1 ? Double.NEGATIVE_INFINITY : motherPost[j] + fatherCum[i + 1];
          nextFather[i]++;
          iterations++;
        }
      }
    }
    final int totalSize = fatherHyp.length * motherHyp.length;
    int fatherRow;

    while ((fatherRow = getNextFatherHyp(fatherPost, motherPost, nextFather)) != -1 && keepGoing(fatherHyp, motherHyp, faHypRem, moHypRem)) {
      final int motherCol = nextFather[fatherRow];
      nextFather[fatherRow]++;

      final int currentFather = fatherHyp[fatherRow];
      final int currentMother = motherHyp[motherCol];

      computeParentMarginal(isDefault, ap, haploidSize, fatherPost[fatherRow], motherPost[motherCol], currentFather, currentMother);

      faHypRem[fatherRow] = motherCol == motherCum.length - 1 ? Double.NEGATIVE_INFINITY : fatherPost[fatherRow] + motherCum[motherCol + 1];
      moHypRem[motherCol] = fatherRow == fatherCum.length - 1 ? Double.NEGATIVE_INFINITY : motherPost[motherCol] + fatherCum[fatherRow + 1];
      iterations++;

    }
//    Diagnostic.developerLog("Broke after " + iterations + "/" + totalSize);
    if (!SQUARE_FIRST) {
      // If we haven't done work for the ref:ref index already do it now.
      if (nextFather[fatherReferenceIndex] <= motherReferenceIndex) {
        computeParentMarginal(isDefault, ap, haploidSize, fatherPost[fatherReferenceIndex], motherPost[motherReferenceIndex], mFatherHypotheses.reference(), mMotherHypotheses.reference());
      }
    }
//    Diagnostic.developerLog("FatherRef: " + fatherReferenceIndex + ", MotherRef: " + motherReferenceIndex + " nextFather: " + Arrays.toString(nextFather));

    HIST.increment(100 * iterations / totalSize);
  }

  private boolean keepGoing(int[] fatherHyp, int[] motherHyp, double[] faHypRem, double[] moHypRem) {
    if (!TERMINATE) {
      return true;
    }
    final int refHyp = mFatherHypotheses.reference();
    assert refHyp == mMotherHypotheses.reference();
    if (refHyp != -1) {
      final double fatherRef = mFatherMarginal[refHyp];
      final double motherRef = mMotherMarginal[refHyp];
      if (fatherRef == Double.NEGATIVE_INFINITY || motherRef == Double.NEGATIVE_INFINITY) {
        return true;
      }
      if (referenceNeedsMoreAccuracy(fatherHyp, faHypRem, fatherRef, refHyp, mFatherMarginal)) {
        return true;
      }
      if (referenceNeedsMoreAccuracy(motherHyp, moHypRem, motherRef, refHyp, mMotherMarginal)) {
        return true;
      }
    }

    if (scoreNeedsMoreAccuracy(fatherHyp, faHypRem, mFatherMarginal)) {
      return true;
    }
    if (scoreNeedsMoreAccuracy(motherHyp, moHypRem, mMotherMarginal)) {
      return true;
    }
    return false;
  }

  /**
   * @param hypothesisLookup array containing the hypothesis corresponding to each entry in the remainders (maps from index in hypothesis remainders to index in marginals)
   * @param hypothesisRemainders un computed portion of each hypothesis
   * @param marginals the marginal array for this individual
   * @return true if more calculation is required for the posterior score
   */
  private boolean scoreNeedsMoreAccuracy(int[] hypothesisLookup, double[] hypothesisRemainders, double[] marginals) {
    double bestFather = Double.NEGATIVE_INFINITY;
    double bestFatherRemainder = Double.NEGATIVE_INFINITY;
    double fatherRemainder = Double.NEGATIVE_INFINITY;
    double fatherSum = Double.NEGATIVE_INFINITY;
    for (int i = 0; i < hypothesisRemainders.length; i++) {
      final double aHypRem = hypothesisRemainders[i];
      final double fatherScore = marginals[hypothesisLookup[i]];
      if (fatherScore > bestFather) {
        fatherRemainder = VariantUtils.logSumApproximation(fatherRemainder, bestFatherRemainder);
        fatherSum = VariantUtils.logSumApproximation(fatherSum, bestFather);
        bestFather = fatherScore;
        bestFatherRemainder = aHypRem;
      } else {
        fatherRemainder = VariantUtils.logSumApproximation(aHypRem, fatherRemainder);
        fatherSum = VariantUtils.logSumApproximation(fatherScore, fatherSum);
      }
      //        remainder = VariantUtils.logSumApproximation(remainder, aHypRem);
    }

    return bestFather == Double.NEGATIVE_INFINITY || (bestFather - bestFatherRemainder) <= TERM_THRESHOLD
      || fatherSum == Double.NEGATIVE_INFINITY || (fatherSum - fatherRemainder) <= TERM_THRESHOLD
      || (bestFather - fatherRemainder) <= TERM_THRESHOLD;
  }

  /**
   * @param hypothesisLookup array containing the hypothesis corresponding to each entry in the remainders (maps from index in hypothesis remainders to index in marginals)
   * @param hypothesisRemainders un computed portion of each hypothesis
   * @param referenceHypothesisScore computed score for the reference hypothesis
   * @param referenceHypothesis reference hypothesis
   * @param marginals the marginal array for this individual
   * @return true if more calculation is required for the non identity score of this individual
   */
  private boolean referenceNeedsMoreAccuracy(int[] hypothesisLookup, double[] hypothesisRemainders, double referenceHypothesisScore, int referenceHypothesis, double[] marginals) {
    double refFatherRemainder = Double.NEGATIVE_INFINITY;
    double fatherRefRemainder = Double.NEGATIVE_INFINITY;
    double fatherRefSum = Double.NEGATIVE_INFINITY;
    for (int i = 0; i < hypothesisRemainders.length; i++) {
      final double aHypRem = hypothesisRemainders[i];
      if (hypothesisLookup[i] == referenceHypothesis) {
        refFatherRemainder = aHypRem;
      } else {
        final double fatherScore = marginals[hypothesisLookup[i]];
        fatherRefRemainder = VariantUtils.logSumApproximation(aHypRem, fatherRefRemainder);
        fatherRefSum = VariantUtils.logSumApproximation(fatherScore, fatherRefSum);
      }
      //        remainder = VariantUtils.logSumApproximation(remainder, aHypRem);
    }
    return (referenceHypothesisScore - refFatherRemainder) <= TERM_THRESHOLD
      || fatherRefSum == Double.NEGATIVE_INFINITY || (fatherRefSum - fatherRefRemainder) <= TERM_THRESHOLD;
  }

  private void computeParentMarginal(boolean isDefault, AlleleProbability ap, int haploidSize, double fatherPost, double motherPost, int currentFather, int currentMother) {
    final double marginal;
    if (isDefault) {
      final double generatedPrior = ap.probabilityLn(mMaximalCode, mReferenceHypothesis, currentFather, currentMother);
      // p prime in theory doc
      final int alleles = Utils.numberAlleles(haploidSize, mMaximalCode, mReferenceHypothesis, currentFather, currentMother);
      marginal =
          -BinomialSpecial.logBinomial(haploidSize, alleles)
              + mParams.getAlleleFrequencyLnProbability(alleles)
              + fatherPost + motherPost
              + generatedPrior;
    } else {
      final double fatherMarginal = currentFather == MISSING_HYPOTHESIS ? 0 : fatherPost;
      final double motherMarginal = currentMother == MISSING_HYPOTHESIS ? 0 : motherPost;
      marginal = fatherMarginal + motherMarginal;
    }

//    System.err.println(currentFather + " " + currentMother + " " + marginal);
    calculateMarginals(marginal, currentFather, currentMother);
  }

  private double[] makeCumulative(double[] fatherPost) {
    final double[] fatherCum = Arrays.copyOf(fatherPost, fatherPost.length);
    for (int i = fatherCum.length - 2; i >= 0; i--) {
      fatherCum[i] = VariantUtils.logSumApproximation(fatherCum[i + 1], fatherCum[i]);
    }
    return fatherCum;
  }

  static final SpyHistogram HIST = new SpyHistogram("SuperAlleles", 100);

  private int getNextFatherHyp(double[] fatherPost, double[] motherPost, int[] next) {
    int best = -1;
    double bestScore = Double.NEGATIVE_INFINITY;
    for (int i = 0; i < fatherPost.length; i++) {
      final int j = next[i];
      if (j >= motherPost.length) {
        continue;
      }
      final double currentScore = fatherPost[i] + motherPost[j];
      if (currentScore > bestScore) {
        bestScore = currentScore;
        best = i;
      }
      // Optional speedup
      /*
      if (next[i] == 0) {
        break;
      }
      */
    }
    return best;
  }

  double calculateMarginals(double marginal, int father, int mother) {
    return calculateMarginalsSlowly(marginal, father, mother, new Stack<Integer>(), 0);
  }

  protected static AlleleProbability getAlleleProbability(Ploidy father, Ploidy mother) {
    if (father == Ploidy.NONE) {
      // Father does not have this sequence
      if (mother == Ploidy.NONE) {
        throw new UnsupportedOperationException("Mother none, father none, not supported.");
      }
      if (mother != Ploidy.HAPLOID) {
        throw new UnsupportedOperationException("Mother diploid, father none, not supported.");
      }
      return AlleleProbabilityHN.SINGLETON_NH;
    }
    if (father == Ploidy.HAPLOID) {
      if (mother == Ploidy.NONE) {
        // Mother does not have this sequence
        return AlleleProbabilityHN.SINGLETON_HN;
      }
      if (mother == Ploidy.HAPLOID) {
        return AlleleProbabilityHH.SINGLETON;
      }
      return AlleleSetProbabilityHaploid.SINGLETON_HD;
    }
    if (mother == Ploidy.NONE) {
      throw new UnsupportedOperationException("Father diploid, mother none, not supported.");
    }
    if (mother == Ploidy.HAPLOID) {
      return AlleleSetProbabilityHaploid.SINGLETON_DH;
    }
    return AlleleSetProbabilityDiploid.SINGLETON;
  }

  void updateGlobalIndentityMarginal(int father, int mother, Stack<Integer> children, double marginal) {
    if (father == MISSING_HYPOTHESIS && mother == MISSING_HYPOTHESIS) {
      throw new RuntimeException("missing sequence from both parents");
    }
    boolean nonIdentity = parentsNonIdentity(father, mother);
    for (int i = 0; !nonIdentity && i < children.size(); i++) {
      final int child = children.get(i);
      if (child != MISSING_HYPOTHESIS && child != mReferenceHypothesis) {
        nonIdentity = true;
        break;
      }
    }
    if (nonIdentity) {
      mNonIdentity = VariantUtils.logSumApproximation(mNonIdentity, marginal);
    } else {
      mIdentity = VariantUtils.logSumApproximation(mIdentity, marginal);
    }
  }

  boolean parentsNonIdentity(int father, int mother) {
    return father != MISSING_HYPOTHESIS && father != mReferenceHypothesis || mother != MISSING_HYPOTHESIS && mother != mReferenceHypothesis;
  }

  private double safeGetCount(final AlleleStatistics<?> counts, final int allele) {
    return allele < counts.getDescription().size() ? counts.count(allele) : 0;
  }

  protected double contraryEvidenceAdjustment(final int fatherHyp, final int motherHyp, final int childHyp) {
    if (mContraryProbabilityLn < 0) { // Efficiency
      // Corresponds to R(H_c | H_m, H_f, E_c, E_m, E_f) in theory document
      final AlleleStatistics<?> fatherCounts = mFather.statistics().counts();
      final AlleleStatistics<?> motherCounts = mMother.statistics().counts();
      final int childA = mMaximalCode.a(childHyp);
      final int childB = mMaximalCode.bc(childHyp);
      final int fatherA = mMaximalCode.a(fatherHyp);
      final int fatherB = mMaximalCode.bc(fatherHyp);
      final int motherA = mMaximalCode.a(motherHyp);
      final int motherB = mMaximalCode.bc(motherHyp);
      if (childHyp != fatherHyp && childHyp != motherHyp) { // No adjustment needed in case where hypotheses are the same
        // This needs to be careful to ensure each piece of contrary evidence is only counted once
        double contraryCount = 0;
        if (childA != fatherA && childA != fatherB && childA != motherA && childA != motherB) {
          // child A allele is not present in hypothesis for either parent
          contraryCount += safeGetCount(fatherCounts, childA) + safeGetCount(motherCounts, childA);
        }
        if (childB != childA && childB != fatherA && childB != fatherB && childB != motherA && childB != motherB) {
          // child heterozygous B allele is not present in hypothesis for either parent
          contraryCount += safeGetCount(fatherCounts, childB) + safeGetCount(motherCounts, childB);
        }
        return mContraryProbabilityLn * contraryCount;
      }
    }
    return 0; // ln(1)
  }

  private double calculateMarginalsSlowly(double marginal, int father, int mother, Stack<Integer> childHyps, int childCount) {
    double result = Double.NEGATIVE_INFINITY;
    if (childCount == mChildren.size()) {
      updateGlobalIndentityMarginal(father, mother, childHyps, marginal);
      //update marginals.
      if (father != MISSING_HYPOTHESIS) {
        mFatherMarginal[father] = VariantUtils.logSumApproximation(mFatherMarginal[father], marginal);
      }
      if (mother != MISSING_HYPOTHESIS) {
        mMotherMarginal[mother] = VariantUtils.logSumApproximation(mMotherMarginal[mother], marginal);
      }
      result = marginal;
      for (int i = 0; i < childHyps.size(); i++) {
        final double[] marginals = mChildMarginal.get(i);
        if (marginals.length > 0) {
          final int cat = childHyps.get(i);
          marginals[cat] = VariantUtils.logSumApproximation(marginals[cat], marginal);
          if (mChildDenovoMarginal != null) {
            if (isDenovo(father, mother, i, cat)) {
              final double denovoMarginal = mChildDenovoMarginal[i];
              mChildDenovoMarginal[i] = VariantUtils.logSumApproximation(denovoMarginal, marginal);
            } else {
              final double nonDenovoMarginal = mChildNonDenovoMarginal[i];
              mChildNonDenovoMarginal[i] = VariantUtils.logSumApproximation(nonDenovoMarginal, marginal);
            }
          }
        }
      }
    } else {
      final ModelInterface<?> child = mChildren.get(childCount);
      final Hypotheses<?> childHypotheses = mHypotheses.get(child);
      final int childHypothesesSize = childHypotheses.size();
      // All the references should be the same, so just use the child's since we have it
      final MendelianAlleleProbability map = MendelianAlleleProbabilityFactory.COMBINED.getMendelianAlleleProbability(mFatherPloidy, mMotherPloidy, childHypotheses.ploidy(), mLogDenovoRefPrior, mLogDenovoNonrefPrior, child.reference());
      if (childHypothesesSize == 0) {
        //near as I can tell this will always be zero currently. MISSING_HYPOTHESIS is ignored
        final double prior = map.probabilityLn(mMaximalCode, father, mother, MISSING_HYPOTHESIS);
        final double childMarginal = marginal + prior;
        //dunno about hear
        childHyps.push(MISSING_HYPOTHESIS);
        result = VariantUtils.logSumApproximation(result, calculateMarginalsSlowly(childMarginal, father, mother, childHyps, childCount + 1));
        childHyps.pop();
      } else {

        for (int i = 0; i < childHypothesesSize; i++) {
          final double prior = map.probabilityLn(mMaximalCode, father, mother, i);
          if (Double.isInfinite(prior)) {
            continue;
          }
          final double childMarginal = marginal + child.posteriorLn0(i) + prior + contraryEvidenceAdjustment(father, mother, i);
          assert !Double.isNaN(childMarginal) : marginal + ":" + child.posteriorLn0(i) + ":" + prior;
          childHyps.push(i);
          result = VariantUtils.logSumApproximation(result, calculateMarginalsSlowly(childMarginal, father, mother, childHyps, childCount + 1));
          childHyps.pop();
        }
      }
    }
    return result;
  }

  protected static void markDenovo(HypothesisScore score, double denovo, double nonDenovo) {
    if (denovo > nonDenovo) {
      score.setDenovo(DeNovoStatus.IS_DE_NOVO);
    }
    score.setDeNovoPosterior(denovo - nonDenovo);
  }

  @Override
  public double getNonIdentityPosterior() {
    return mNonIdentity - mIdentity;
  }

  final boolean isDenovo(int father, int mother, int childId, int childHypothesis) {
    final ModelInterface<?> childModel = mChildren.get(childId);
    final Hypotheses<?> hypotheses = mHypotheses.get(childModel);
    final MendelianAlleleProbability al = MendelianAlleleProbabilityFactory.COMBINED.getMendelianAlleleProbability(mFatherPloidy, mMotherPloidy, hypotheses.ploidy(), mLogDenovoRefPrior, mLogDenovoNonrefPrior, hypotheses.reference());
    return al.isDenovo(mMaximalCode, father, mother, childHypothesis);
  }
}
