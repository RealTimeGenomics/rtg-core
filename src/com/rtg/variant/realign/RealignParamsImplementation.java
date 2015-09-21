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

package com.rtg.variant.realign;

import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;
import com.rtg.variant.AbstractMachineErrorParams;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;

/**
 */
public final class RealignParamsImplementation extends IntegralAbstract implements RealignParams {

  private static final double DEFAULT_DECAY = Math.log(0.2);

  /** Key for asking about the CG overlap. */
  public static final int CG_OVERLAP = 0;

  /** Key for asking about the CG small gap. */
  public static final int CG_SMALL_GAP = 1;

  /** Key for asking about the CG large gap. */
  public static final int CG_LARGE_GAP = 2;

  /** Key for asking about the CG V2 overlap. */
  public static final int CG_OVERLAP2 = 3;

  /**
   * Compute a heuristic decay constant for the distribution.
   * @param distr probability distribution with entry 0 equal to 0.0
   * @return natural log of a decay constant.
   */
  static double weightedDecay(final double[] distr) {
    //should be robust against almost any distribtuion
    if (distr.length < 3) {
      return DEFAULT_DECAY;
    }
    double wtdSum = 0.0;
    double sum = 0.0;
    for (int i = 2; i < distr.length; i++) {
      final double di = distr[i];
      final double dim = distr[i - 1];
      assert dim >= 0.0 && dim <= 1.0;
      if (dim > 0 && di > 0) {
        final double ratio = di / dim;
        final double d = dim * Math.log(ratio);
        //System.err.println("i=" + i + " di=" + dim + " ratio=" + ratio + " d=" + d);
        wtdSum += d;
        sum += dim;
      }
    }
    if (sum == 0.0) {
      return DEFAULT_DECAY;
    }
    assert sum > 0.0;
    final double res = wtdSum / sum;
    assert !Double.isInfinite(res) && !Double.isNaN(res);
    return res;
  }

  /**
   * Compute the mean value of the distribution.
   * @param distr whose mean is to be computed.
   * @return mean.
   */
  static double mean(final double[] distr) {
    double sum = 0.0;
    for (int i = 1; i < distr.length; i++) {
      sum += i * distr[i];
    }
    assert 1.0 <= sum && sum <= (distr.length - 1) && !Double.isNaN(sum) : sum;
    return sum;
  }

  private final double mDeleteOpenLn;

  private final double mDeleteExtendLn;

  private final double mInsertOpenLn;

  private final double mInsertExtendLn;

  private final double mMatchLn;

  private final double mMisMatchLn;

  private final boolean mCG;

  private final int[] mGapStart;
  private final double[][] mGapDistributions;

  private void check(final double x) {
    Exam.assertTrue(x <= 0.0 && !Double.isNaN(x));
  }

  @Override
  public boolean integrity() {
    check(mDeleteOpenLn);
    check(mDeleteExtendLn);
    check(mInsertOpenLn);
    check(mInsertExtendLn);
    check(mMatchLn);
    check(mMisMatchLn);
    Exam.assertEquals(1.0, Math.exp(mMatchLn) + Math.exp(mMisMatchLn), 0.000001);
    return true;
  }


  /**
   * @param params Priors from a priors file.
   */
  public RealignParamsImplementation(final AbstractMachineErrorParams params) {
    mDeleteOpenLn = Math.log(params.errorDelEventRate());
    mDeleteExtendLn = weightedDecay(params.errorDelDistribution());
    mInsertOpenLn = Math.log(params.errorInsEventRate());
    mInsertExtendLn = weightedDecay(params.errorInsDistribution());
    final double misMatch = mean(params.errorMnpDistribution()) * params.errorMnpEventRate();
    mMisMatchLn = Math.log(misMatch);
    mMatchLn = Math.log(1.0 - misMatch);
    mCG = params.isCG();
    mGapDistributions = new double[4][];
    mGapStart = new int[] {-4, 0, 4, -7};
    logify(CG_OVERLAP, params.overlapDistribution());
    logify(CG_SMALL_GAP, params.smallGapDistribution());
    logify(CG_LARGE_GAP, params.gapDistribution());
    logify(CG_OVERLAP2, params.overlapDistribution2());
  }

  private void logify(final int whichGap, final double[] probs) {
    // we strip leading and trailing zeroes
    //System.err.println("logify " + whichGap + " in=" + Arrays.toString(probs));
    int start = 0;
    while (start < probs.length && probs[start] == 0.0) {
      start++;
      mGapStart[whichGap]++;
    }
    int end = probs.length - 1;
    while (end > start && probs[end] == 0.0) {
      end--;
    }
    mGapDistributions[whichGap] = new double[end + 1 - start];
    for (int i = start; i <= end; i++) {
      mGapDistributions[whichGap][i - start] = Math.log(probs[i]);
    }
    //System.err.println("logify " + whichGap + " out=" + Arrays.toString(mGapDistributions[whichGap]));
  }

  @Override
  public double[][] gapDistributionPoss(final PossibilityArithmetic arith) {
    return gapDistributionPoss(mGapDistributions, arith);
  }

  /**
   * @param gapDistributions gap distributions in log format.
   * @param arith arithmetic delegate.
   * @return an array of possibilities (as determined by arith) for each CG gap.
   */
  public static double[][] gapDistributionPoss(final double[][] gapDistributions, final PossibilityArithmetic arith) {
    final double[][] gdp = new double[gapDistributions.length][];
    for (int i = 0; i < gapDistributions.length; i++) {
      gdp[i] = new double[gapDistributions[i].length];
      for (int j = 0; j < gapDistributions[i].length; j++) {
        gdp[i][j] = arith.ln2Poss(gapDistributions[i][j]);
      }
    }
    return gdp;
  }

  @Override
  public double deleteExtendLn() {
    return mDeleteExtendLn;
  }

  @Override
  public double deleteOpenLn() {
    return mDeleteOpenLn;
  }

  @Override
  public double insertExtendLn() {
    return mInsertExtendLn;
  }

  @Override
  public double insertOpenLn() {
    return mInsertOpenLn;
  }

  @Override
  public double matchLn() {
    return mMatchLn;
  }

  @Override
  public double misMatchLn() {
    return mMisMatchLn;
  }

  @Override
  public boolean completeGenomics() {
    return mCG;
  }

  @Override
  public int gapStart(final int gap) {
    return mGapStart[gap];
  }

  @Override
  public int gapEnd(final int gap) {
    return mGapStart[gap] + mGapDistributions[gap].length - 1;
  }

  @Override
  public double gapFreqLn(final int gap, final int width) {
    return mGapDistributions[gap][width - gapStart(gap)];
  }

  @Override
  public void toString(final StringBuilder sb) {
    sb.append(false); //TODO
  }

}
