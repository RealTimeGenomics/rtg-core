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

import java.util.List;

import com.rtg.relation.Family;
import com.rtg.variant.GenomePriorParams;
import com.rtg.variant.bayes.Hypotheses;
import com.rtg.variant.bayes.ModelInterface;
import com.rtg.variant.bayes.multisample.Utils;
import com.rtg.variant.util.VariantUtils;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;

/**
 */
public class FastDiseasedFamilyPosterior extends DiseasedFamilyPosterior {
  static final PossibilityArithmetic SINGLETON = com.rtg.variant.util.arithmetic.LogApproximatePossibility.SINGLETON;
  static final BitSet BIT_SET = BitSet.DNA_SET;

  private final int[] mCodeToBits;
  private final boolean[] mChildDiseased;

  //TODO add flag for recessive vs dominant
  FastDiseasedFamilyPosterior(final GenomePriorParams params, final Family family, Hypotheses<?> commonHypotheses, HypothesesDisease diseaseHypotheses, final List<ModelInterface<?>> models) {
    super(params, family, commonHypotheses, diseaseHypotheses, models);
    mCodeToBits = new int[mChildren.get(0).size()];
    for (int k = 0; k < mCodeToBits.length; k++) {
      mCodeToBits[k] = code2Bits(k);
    }
    mChildDiseased = new boolean[mChildren.size()];
    for (int k = 0; k < mChildDiseased.length; k++) {
      mChildDiseased[k] = mFamily.isDiseased(Family.FIRST_CHILD_INDEX + k);
    }
  }

  private int code2Bits(final int code) {
    if (mCode.homozygous(code)) {
      return BIT_SET.toSet(mCode.a(code));
    } else {
      return BIT_SET.toSet(mCode.a(code), mCode.b(code));
    }
  }


  @Override
  void calculateMarginals(final double b, final int father, final int mother) {
    final int children = mChildren.size();
    final WeightedLattice[] weightedChildren = new WeightedLattice[children];
    for (int i = 0; i < children; i++) {
      weightedChildren[i] = new DefaultWeightedLattice(SINGLETON, BIT_SET);
      final ModelInterface<?> child = mChildren.get(i);
      final boolean diseased = mChildDiseased[i];
      for (int j = 0; j < child.size(); j++) {
        final double mendelian = MendelianAlleleProbabilityDiploid.SINGLETON.probabilityLn(mCode, father, mother, j);
        final int kidSet = mCodeToBits[j];
        final int set = diseased ? kidSet : BIT_SET.complement(kidSet);
        weightedChildren[i].set(set, child.posteriorLn0(j) + mendelian);
      }
    }

    final int diseasedParentSet;
    final int healthyParentSet;
    if (mFamily.isDiseased(Family.FATHER_INDEX)) {
      assert !mFamily.isDiseased(Family.MOTHER_INDEX);
      diseasedParentSet = code2Bits(father);
      healthyParentSet = code2Bits(mother);
    } else {
      assert mFamily.isDiseased(Family.MOTHER_INDEX);
      diseasedParentSet = code2Bits(mother);
      healthyParentSet = code2Bits(father);
    }

    final int initSet = diseasedParentSet & BIT_SET.complement(healthyParentSet) & BIT_SET.complement(BIT_SET.toSet(mHypotheses.reference()));

    final WeightedLattice init = new DefaultWeightedLattice(SINGLETON, BIT_SET);
    init.set(initSet, SINGLETON.one());
    final SFunction sf = new SFunction(init, weightedChildren); //TODO use version that does recessive rather than dominant

    // Compute C(H_f,H_m,H_d) for each H_d, pos 0 is NO_DISEASE
    final double[] c = new double[mDiseaseHypotheses.size()];
    for (int d = 0; d < mDiseaseHypotheses.size(); d++) {
      final int alleleCount = Utils.numberAllelesExclude(mHypotheses.description().size(), mCode, d - 1, father, mother, mHypotheses.reference());
      final double pPrime = mParams.getAlleleFrequencyLnProbability(alleleCount);
      c[d] = pPrime + mDiseaseHypotheses.prior(d);
    }

    final WeightedLattice sfAll = sf.all();
    final WeightedLattice.Visitor vall = new WeightedLattice.Visitor() {
      @Override
      public void value(int set, double weight) {
        final double v = SINGLETON.poss2Ln(weight);
        if (set == 0) {
          final double t = v + c[0] + b;
          mDiseaseMarginal[0] = VariantUtils.logSumApproximation(mDiseaseMarginal[0], t);
          mMotherMarginal[mother][0] = VariantUtils.logSumApproximation(mMotherMarginal[mother][0], t);
          mFatherMarginal[father][0] = VariantUtils.logSumApproximation(mFatherMarginal[father][0], t);
        } else {
          final BitSet.Visitor vbit = new BitSet.Visitor() {
            @Override
            public void visit(final int x) {
              final int d = x + 1;
              final double t = v + c[d] + b;
              mDiseaseMarginal[d] = VariantUtils.logSumApproximation(mDiseaseMarginal[d], t);
              mMotherMarginal[mother][d] = VariantUtils.logSumApproximation(mMotherMarginal[mother][d], t);
              mFatherMarginal[father][d] = VariantUtils.logSumApproximation(mFatherMarginal[father][d], t);
            }
          };
          BIT_SET.visit(vbit, set);
        }
      }
    };
    sfAll.visit(vall);

    //find best disease
    int be = 0;
    double bestWeight = mDiseaseMarginal[0];
    for (int i = 1; i < mDiseaseMarginal.length; i++) {
      if (mDiseaseMarginal[i] > bestWeight) {
        be = i;
        bestWeight = mDiseaseMarginal[i];
      }
    }
    final int bestDisease = be;
    final double cd = b + c[bestDisease];
    final double c0 = b + c[0];
    final int bestSet = bestDisease == 0 ? 0 : BIT_SET.toSet(bestDisease - 1);

    for (int ch = 0; ch < children; ch++) {
      final double[][] childMarginal = mChildMarginal.get(ch);
      final WeightedLattice sfChild = sf.excludeChild(ch);
      final WeightedLattice child = weightedChildren[ch];
      final boolean diseased = mChildDiseased[ch];

      final ModelInterface<?> kittens = mChildren.get(ch); //TODO can this be optimized?
      for (int i = 0; i < kittens.size(); i++) {
        final int hyp = i; //We are speechless
        final int kitSet0 = code2Bits(i);
        final int kitSet = diseased ? kitSet0 : BIT_SET.complement(kitSet0);
        final double w = SINGLETON.poss2Ln(child.get(kitSet));
        final WeightedLattice.Visitor cv = new WeightedLattice.Visitor() {
          @Override
          public void value(final int set, final double weight) {
            final int s = set & kitSet;
            if (s == 0 && bestDisease == 0) {
              final double v = SINGLETON.poss2Ln(weight);
              final double t = c0 + v + w;
              childMarginal[hyp][0] = VariantUtils.logSumApproximation(childMarginal[hyp][0], t);
            } else if ((s & bestSet) != 0) {
              final double v = SINGLETON.poss2Ln(weight);
              final double t = cd + v + w;
              childMarginal[hyp][bestDisease] = VariantUtils.logSumApproximation(childMarginal[hyp][bestDisease], t);
            }
          }
        };
        sfChild.visit(cv);
      }
    }
  }
}
