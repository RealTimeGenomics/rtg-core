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

package com.rtg.variant.bayes.multisample.cancer;

import static com.rtg.util.StringUtils.LS;

import java.util.Arrays;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.util.StringUtils;
import com.rtg.util.Utils;
import com.rtg.util.format.FormatReal;
import com.rtg.variant.bayes.ArrayGenotypeMeasure;
import com.rtg.variant.bayes.GenotypeMeasure;
import com.rtg.variant.bayes.Hypotheses;
import com.rtg.variant.bayes.NoNonIdentityMeasure;
import com.rtg.variant.bayes.Statistics;
import com.rtg.variant.util.arithmetic.LogApproximatePossibility;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;

/**
 * Common part of posterior calculations from pure and contaminated somatic calling.
 */
@TestClass(value = "com.rtg.variant.bayes.multisample.cancer.SomaticPosteriorPureTest")
public abstract class AbstractSomaticPosterior {

  protected final PossibilityArithmetic mArithmetic = LogApproximatePossibility.SINGLETON;
  protected final Hypotheses<?> mNormalHypotheses;
  protected final Hypotheses<?> mCancerHypotheses;
  protected double mEqual = Double.NEGATIVE_INFINITY;
  protected double mNotEqual = Double.NEGATIVE_INFINITY;
  protected final double[][] mPosterior;
  protected final double[] mNormalMarginal;
  protected final double[] mCancerMarginal;
  protected int mBestNormal;
  protected int mBestCancer;
  private final double mPhi;
  private final double mPsi;

  /**
   * @param normalHypotheses for normal sample hypotheses
   * @param cancerHypotheses the cancer sample hypotheses
   * @param phi probability of seeing contrary evidence in the original
   * @param psi probability of seeing contrary evidence in the derived
   */
  public AbstractSomaticPosterior(final Hypotheses<?> normalHypotheses, final Hypotheses<?> cancerHypotheses, final double phi, final double psi) {
    mNormalHypotheses = normalHypotheses;
    mCancerHypotheses = cancerHypotheses;
    mPosterior = new double[normalHypotheses.size()][cancerHypotheses.size()];
    mNormalMarginal = new double[normalHypotheses.size()];
    mCancerMarginal = new double[cancerHypotheses.size()];
    Arrays.fill(mNormalMarginal, Double.NEGATIVE_INFINITY);
    Arrays.fill(mCancerMarginal, Double.NEGATIVE_INFINITY);
    mPhi = mArithmetic.prob2Poss(phi);
    mPsi = mArithmetic.prob2Poss(psi);
  }

  protected double logSum(final double x, final double y) {
    return mArithmetic.add(x, y);
    //return VariantUtils.logSum(x, y);  // more accurate
  }

  protected boolean isSame(final int normalHyp, final int cancerHyp) {
    return normalHyp == cancerHyp;
  }

  /**
   * Finish the construction after values put into <code>mPosterior</code>
   */
  protected void postConstruction() {
    // Calculation marginals and equal/not-equal probabilities
    for (int normal = 0; normal < mNormalHypotheses.size(); ++normal) {
      assert mPosterior[normal].length == mCancerHypotheses.size();
      for (int cancer = 0; cancer < mCancerHypotheses.size(); ++cancer) {
        final double p = mPosterior[normal][cancer];
        if (isSame(normal, cancer)) {
          mEqual = logSum(mEqual, p);
        } else {
          mNotEqual = logSum(mNotEqual, p);
        }
        mNormalMarginal[normal] = logSum(mNormalMarginal[normal], p);
        mCancerMarginal[cancer] = logSum(mCancerMarginal[cancer], p);
      }
    }

    mBestNormal = findBest(mNormalMarginal);
    mBestCancer = findBest(mCancerMarginal);
  }

  protected void adjustPosterior(final int normal, final int cancer, final double contraryProb, final double count) {
    mPosterior[normal][cancer] = mArithmetic.multiply(mPosterior[normal][cancer], mArithmetic.pow(contraryProb, count));
  }

  protected int normalHypToAlleleBits(final int normalHyp) {
    final int normalA = mNormalHypotheses.code().a(normalHyp);
    final int normalB = mNormalHypotheses.code().bc(normalHyp);
    return (1 << normalA) | (1 << normalB);
  }

  protected int cancerHypToAlleleBits(final int cancerHyp) {
    return normalHypToAlleleBits(cancerHyp); // By default they have the same hypotheses as the normal
  }

  protected final void contraryEvidenceAdjustment(final Statistics<?> normalStats, final Statistics<?> cancerStats) {
    // Corresponds to R(H_c | H_n, E_c, E_n) in theory document
    for (int normal = 0; normal < mNormalHypotheses.size(); ++normal) {
      final int normalAlleles = normalHypToAlleleBits(normal);
      for (int cancer = 0; cancer < mCancerHypotheses.size(); ++cancer) {
        final int cancerAlleles = cancerHypToAlleleBits(cancer);
        for (int alleleBit = 1, allele = 0; allele < mNormalHypotheses.description().size(); ++allele, alleleBit <<= 1) {
          if ((cancerAlleles & alleleBit) != 0) {
            if ((normalAlleles & alleleBit) == 0) {
              // allele present in cancer but not normal, any normal evidence is contrary
              adjustPosterior(normal, cancer, mPhi, normalStats.counts().count(allele));
            }
          } else {
            if ((normalAlleles & alleleBit) != 0) {
              // allele present in normal but not in cancer, any cancer evidence is contrary
              adjustPosterior(normal, cancer, mPsi, cancerStats.counts().count(allele));
            }
          }
        }
      }
    }
  }

  int findBest(double[] marginals) {
    double bestScore = mArithmetic.zero();
    int besti = -1;
    for (int i = 0; i < marginals.length; ++i) {
      final double p = marginals[i];
      if (mArithmetic.gt(p, bestScore)) {
        bestScore = p;
        besti = i;
      }
    }
    return besti;
  }

  @Override
  public String toString() {
    final StringBuilder sb = new StringBuilder();
    final FormatReal fmt = new FormatReal(4, 3);
    final int pad = mNormalHypotheses.maxNameLength();
    final int col = Math.max(8, mCancerHypotheses.maxNameLength() + 1);
    sb.append(StringUtils.padLeft("", pad));
    for (int cancerHyp = 0; cancerHyp < mCancerHypotheses.size(); ++cancerHyp) {
      sb.append(StringUtils.padLeft(mCancerHypotheses.name(cancerHyp), col));
    }
    sb.append(LS);
    for (int normalHyp = 0; normalHyp < mNormalHypotheses.size(); ++normalHyp) {
      sb.append(StringUtils.padLeft(mNormalHypotheses.name(normalHyp), pad));
      for (int cancerHyp = 0; cancerHyp < mCancerHypotheses.size(); ++cancerHyp) {
        sb.append(StringUtils.padLeft(fmt.format(mPosterior[normalHyp][cancerHyp]), col));
      }
      sb.append(fmt.format(mNormalMarginal[normalHyp])).append(LS);
    }
    sb.append(StringUtils.padLeft("", pad));
    for (int cancerHyp = 0; cancerHyp < mCancerHypotheses.size(); ++cancerHyp) {
      sb.append(StringUtils.padLeft(fmt.format(mCancerMarginal[cancerHyp]), col));
    }
    sb.append(LS);
    sb.append("best[").append(mBestNormal).append(",").append(mBestCancer).append("]=").append(Utils.realFormat(mPosterior[mBestNormal][mBestCancer], 3)).append(LS);
    sb.append("equal=").append(Utils.realFormat(mEqual, 3)).append("  notequal=").append(Utils.realFormat(mNotEqual, 3)).append(LS);
    return sb.toString();
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
    for (int normalHyp = 0; normalHyp < mPosterior.length; ++normalHyp) {
      for (int cancerHyp = 0; cancerHyp < mPosterior[normalHyp].length; ++cancerHyp) {
        if (normalHyp != mBestNormal || cancerHyp != mBestCancer) {
          others = logSum(others, mPosterior[normalHyp][cancerHyp]);
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
    return new NoNonIdentityMeasure(new ArrayGenotypeMeasure(mArithmetic, mNormalMarginal, mNormalHypotheses));
  }

  /**
   * @return <code>ln(P / (1-P))</code> where P is the posterior for the most likely cancer call.
   */
  protected GenotypeMeasure cancerMeasure() {
    return new NoNonIdentityMeasure(new ArrayGenotypeMeasure(mArithmetic, mCancerMarginal, mCancerHypotheses));
  }

}
