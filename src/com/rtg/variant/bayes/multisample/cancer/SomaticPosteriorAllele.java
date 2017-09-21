/*
 * Copyright (c) 2017. Real Time Genomics Limited.
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

import com.rtg.variant.bayes.Hypotheses;
import com.rtg.variant.bayes.HypothesesPowerSet;
import com.rtg.variant.bayes.Statistics;

/**
 * Posterior calculations for allele based cancer calling.
 */
public class SomaticPosteriorAllele extends AbstractSomaticPosterior {

  /**
   * @param normalHypotheses hypotheses for normal model
   * @param cancerHypotheses hypotheses for cancer model
   * @param phi probability of seeing contrary evidence in the original
   * @param psi probability of seeing contrary evidence in the derived
   */
  public SomaticPosteriorAllele(final Hypotheses<?> normalHypotheses, final HypothesesPowerSet<?> cancerHypotheses, final double phi, final double psi) {
    super(normalHypotheses, cancerHypotheses, phi, psi);
  }

  @Override
  protected void postConstruction() {
    assert mPosterior.length == mNormalHypotheses.size();
    for (int normal = 0; normal < mPosterior.length; ++normal) {
      for (int cancer = 0; cancer < mPosterior[normal].length; ++cancer) {
        final double p = mPosterior[normal][cancer];
        if ((1 << normal) == cancer) {
          mEqual = logSum(mEqual, p);
        } else {
          mNotEqual = logSum(mNotEqual, p);
        }
        mNormalMarginal[normal] = logSum(mNormalMarginal[normal], p);
        mCancerMarginal[cancer] = logSum(mCancerMarginal[cancer], p);
      }
    }
    mBestNormal = AbstractSomaticPosterior.findBest(mNormalMarginal);
    mBestCancer = AbstractSomaticPosterior.findBest(mCancerMarginal);
  }

  @Override
  protected void contraryEvidenceAdjustment(final Statistics<?> normalStats, final Statistics<?> cancerStats) {
    // todo, parent implementation is wrong for these cancer hyp
    // todo for the moment this is making no adjustment
//    // Corresponds to R(H_c | H_n, E_c, E_n) in theory document
//    for (int normal = 0; normal < mLength; ++normal) {
//      final int normalA = mNormalHypotheses.code().a(normal);
//      final int normalB = mNormalHypotheses.code().bc(normal);
//      for (int cancer = 0; cancer < mLength; ++cancer) {
//        if (normal != cancer) { // No adjustment needed in case where hypotheses are the same
//          final int cancerA = mNormalHypotheses.code().a(cancer);
//          final int cancerB = mNormalHypotheses.code().bc(cancer);
//          double contraryNormalCount = 0;
//          double contraryCancerCount = 0;
//          if (cancerA != normalA && cancerA != normalB) {
//            contraryNormalCount += normalStats.counts().count(cancerA);
//          }
//          if (cancerB != cancerA && cancerB != normalA && cancerB != normalB) { // Only for heterozygous
//            contraryNormalCount += normalStats.counts().count(cancerB);
//          }
//          if (normalA != cancerA && normalA != cancerB) {
//            contraryCancerCount += cancerStats.counts().count(normalA);
//          }
//          if (normalB != normalA && normalB != cancerA && normalB != cancerB) { // Only for heterozygous
//            contraryCancerCount += cancerStats.counts().count(normalB);
//          }
//          mPosterior[normal][cancer] += mPhi * contraryNormalCount + mPsi * contraryCancerCount;
//        }
//      }
//    }
  }
}
