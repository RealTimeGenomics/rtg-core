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

package com.rtg.variant.bayes.multisample.cancer;

import java.util.Arrays;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.util.StringUtils;
import com.rtg.util.Utils;
import com.rtg.util.format.FormatReal;
import com.rtg.util.integrity.IntegralAbstract;
import com.rtg.variant.bayes.ArrayGenotypeMeasure;
import com.rtg.variant.bayes.GenotypeMeasure;
import com.rtg.variant.bayes.Hypotheses;
import com.rtg.variant.bayes.NoNonIdentityMeasure;
import com.rtg.variant.util.arithmetic.LogApproximatePossibility;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;

/**
 */
@TestClass(value = "com.rtg.variant.bayes.multisample.cancer.PosteriorPureTest")
public abstract class AbstractPosterior extends IntegralAbstract {

  PossibilityArithmetic mArithmetic = LogApproximatePossibility.SINGLETON;
  protected double logSum(final double x, final double y) {
    return mArithmetic.add(x, y);
    //return VariantUtils.logSum(x, y);  // more accurate
  }

  protected final Hypotheses<?> mHypotheses;
  protected final int mLength;
  protected double mEqual = Double.NEGATIVE_INFINITY;
  protected double mNotEqual = Double.NEGATIVE_INFINITY;
  protected final double[][] mPosterior;
  protected final double[] mNormalMarginal;
  protected final double[] mCancerMarginal;
  protected int mBestNormal;
  protected int mBestCancer;

  /**
   * @param hypotheses for normal model (cancer is cross-product of this).
   */
  public AbstractPosterior(final Hypotheses<?> hypotheses) {
    mHypotheses = hypotheses;
    mLength = hypotheses.size();
    mPosterior = new double[mLength][mLength];
    mNormalMarginal = new double[mLength];
    mCancerMarginal = new double[mLength];
    Arrays.fill(mNormalMarginal, Double.NEGATIVE_INFINITY);
    Arrays.fill(mCancerMarginal, Double.NEGATIVE_INFINITY);
  }

  /**
   * Finish the construction after values put into <code>mPosterior</code>
   */
  protected final void postConstruction() {
    // now calculate row/column/diagonal/non-diagonal sums.
    for (int i = 0; i < mLength; i++) {
      assert mPosterior[i].length == mLength;
      for (int j = 0; j < mLength; j++) {
        final double p = mPosterior[i][j];
        if (i == j) {
          mEqual = logSum(mEqual, p);
        } else {
          mNotEqual = logSum(mNotEqual, p);
        }
        mNormalMarginal[i] = logSum(mNormalMarginal[i], p);
        mCancerMarginal[j] = logSum(mCancerMarginal[j], p);
      }
    }

    mBestNormal = findBest(mNormalMarginal);
    mBestCancer = findBest(mCancerMarginal);
  }

  private int findBest(double[] marginals) {
    double bestScore = Double.NEGATIVE_INFINITY;
    int besti = -1;
    for (int i = 0; i < mLength; i++) {
      final double p = marginals[i];
      if (p > bestScore) {
        bestScore = p;
        besti = i;
      }
    }
    return besti;
  }

  @Override
  public void toString(StringBuilder sb) {
    final FormatReal fmt = new FormatReal(4, 3);
    final int pad = mHypotheses.nameLength();
    for (int i = 0; i < mLength; i++) {
      sb.append(StringUtils.padLeft(mHypotheses.name(i), pad));
      for (int j = 0; j < mLength; j++) {
        sb.append(fmt.format(mPosterior[i][j]));
      }
      sb.append(fmt.format(mNormalMarginal[i])).append(LS);
    }
    sb.append(StringUtils.padLeft("", pad));
    for (int j = 0; j < mLength; j++) {
      sb.append(fmt.format(mCancerMarginal[j]));
    }
    sb.append(LS);
    sb.append("best[").append(mBestNormal).append(",").append(mBestCancer).append("]=").append(Utils.realFormat(mPosterior[mBestNormal][mBestCancer], 3)).append(LS);
    sb.append("equal=").append(Utils.realFormat(mEqual, 3)).append("  notequal=").append(Utils.realFormat(mNotEqual, 3)).append(LS);
  }

  /**
   * @return position of the most likely hypothesis for the normal genome,
   *     after calculating the joint posteriors.
   */
  protected int bestNormal() {
    return mBestNormal;
  }

  /**
   * @return position of the most likely hypothesis for the cancer genome,
   *     after calculating the joint posteriors.
   */
  protected int bestCancer() {
    return mBestCancer;
  }

  /**
   * @return <code>ln(P/(1-P))</code>, where P is the posterior of the most likely normal and cancer calls.
   */
  protected double posteriorScore() {
    final double best = mPosterior[mBestNormal][mBestCancer];
    double others = Double.NEGATIVE_INFINITY;
    for (int i = 0; i < mLength; i++) {
      for (int j = 0; j < mLength; j++) {
        if (i != mBestNormal || j != mBestCancer) {
          final double pij = mPosterior[i][j];
          others = logSum(others, pij);
        }
      }
    }
    return best - others;
  }

  /**
   * Calculate the log of <code>Neq / (1 - Neq)</code>,
   * which is the same as <code>Neq / Eq</code>.
   * @return the log posterior for the cancer call being different from the normal call.
   */
  protected double ncScore() {
    return mNotEqual - mEqual;
  }

  /**
   * @return true if the cancer call seems to be the same as the normal call.
   */
  protected boolean isSameCall() {
    return mEqual >= mNotEqual;
  }

  /**
   * @return <code>ln(P / (1-P))</code> where P is the posterior for the most likely normal call.
   */
  protected GenotypeMeasure normalMeasure() {
    return new NoNonIdentityMeasure(new ArrayGenotypeMeasure(mArithmetic, mNormalMarginal, mHypotheses));
  }

  /**
   * @return <code>ln(P / (1-P))</code> where P is the posterior for the most likely cancer call.
   */
  protected GenotypeMeasure cancerMeasure() {
    return new NoNonIdentityMeasure(new ArrayGenotypeMeasure(mArithmetic, mCancerMarginal, mHypotheses));
  }

  @Override
  public boolean integrity() {
    return true;
  }

}
