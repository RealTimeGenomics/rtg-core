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

import java.util.Arrays;
import java.util.List;

import com.rtg.relation.Family;
import com.rtg.variant.GenomePriorParams;
import com.rtg.variant.bayes.Hypotheses;
import com.rtg.variant.bayes.ModelInterface;
import com.rtg.variant.bayes.multisample.HaploidDiploidHypotheses;
import com.rtg.variant.util.VariantUtils;

/**
 * Posterior calculation for the family caller.
 */
public final class FastFamilyPosterior extends FamilyPosterior {

  FastFamilyPosterior(Family family, GenomePriorParams params, List<ModelInterface<?>> models, HaploidDiploidHypotheses<?> hypotheses) {
    super(family, params, models, hypotheses);
  }

  @Override
  double calculateMarginals(double marginal, int father, int mother) {
    final int children = mChildren.size();
    final double[] r = new double[children];
    final double[][] rh = new double[children][];

    for (int i = 0; i < children; i++) {
      final ModelInterface<?> child = mChildren.get(i);
      final Hypotheses<?> childHyp = mHypotheses.get(child);
      // All the references should be the same, so just use the child's since we have it
      final MendelianAlleleProbability map = MendelianAlleleProbabilityFactory.COMBINED.getMendelianAlleleProbability(mFatherPloidy, mMotherPloidy, childHyp.ploidy(), mLogDenovoRefPrior, mLogDenovoNonrefPrior, childHyp.reference());
      rh[i] = new double[mChildren.get(i).size()];
      Arrays.fill(rh[i], Double.NEGATIVE_INFINITY);
      double ri = Double.NEGATIVE_INFINITY;
      for (int j = 0; j < rh[i].length; j++) {
        final double mendelian = map.probabilityLn(mMaximalCode, father, mother, j);
        rh[i][j] = child.posteriorLn0(j) + mendelian + contraryEvidenceAdjustment(child, father, mother, j);
        ri = VariantUtils.logSumApproximation(ri, rh[i][j]);
      }
      r[i] = ri;
    }

    final double[] rForward = new double[children + 1];
    final double[] rReverse = new double[children + 1];
    for (int i = 0; i < children; i++) {
      if (mChildren.get(i).size() > 0) {
        rForward[i + 1] = rForward[i] + r[i];
      } else {
        rForward[i + 1] = rForward[i];
      }
    }

    for (int i = 0; i < children; i++) {
      if (mChildren.get(children - i - 1).size() > 0) {
        rReverse[children - i - 1] = rReverse[children - i] + r[children - i - 1];
      } else {
        rReverse[children - i - 1] = rReverse[children - i];
      }
    }
    final double parentMarginal = rForward[children] + marginal;
    if (father != MISSING_HYPOTHESIS) {
      mFatherMarginal[father] = VariantUtils.logSumApproximation(mFatherMarginal[father], parentMarginal);
    }
    if (mother != MISSING_HYPOTHESIS) {
      mMotherMarginal[mother] = VariantUtils.logSumApproximation(mMotherMarginal[mother], parentMarginal);
    }
    if (!parentsNonIdentity(father, mother)) {
      // parents are both ref so compute the probability of all children being ref and add it to the identity marginal
      double ref = marginal;
      for (int i = 0; i < children; i++) {
        if (rh[i].length > 0) {
          ref += rh[i][mChildren.get(i).reference()];
        }
      }

      // the non identity is the total excluding all ref
      final double nonIdentity;
      if (parentMarginal > ref) {
        nonIdentity = VariantUtils.logSubtract(parentMarginal, ref);
      } else {
        nonIdentity = Double.NEGATIVE_INFINITY;
      }

      mIdentity = VariantUtils.logSumApproximation(mIdentity, ref);
      mNonIdentity = VariantUtils.logSumApproximation(mNonIdentity, nonIdentity);
    } else {
      // parents are non identity so the whole marginal is added to non identity marginal
      mNonIdentity = VariantUtils.logSumApproximation(mNonIdentity, parentMarginal);
    }

    for (int child = 0; child < children; child++) {
      final double[] childMarginals = mChildMarginal.get(child);
      for (int i = 0; i < childMarginals.length; i++) {
        final double hypMarginal = marginal + rh[child][i] + rForward[child] + rReverse[child + 1];
//        System.err.println("child: " + child + ", i: " + i + ", hypMarginal: " + hypMarginal);
        childMarginals[i] = VariantUtils.logSumApproximation(childMarginals[i], hypMarginal);
        if (mChildDenovoMarginal != null) {
          if (isDenovo(father, mother, child, i)) {
            mChildDenovoMarginal[child] = VariantUtils.logSumApproximation(mChildDenovoMarginal[child], hypMarginal);
          } else {
            mChildNonDenovoMarginal[child] = VariantUtils.logSumApproximation(mChildNonDenovoMarginal[child], hypMarginal);
          }
        }

      }
    }
    return parentMarginal;
  }

}
