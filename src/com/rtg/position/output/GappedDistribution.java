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
package com.rtg.position.output;

import static com.rtg.position.output.Offset.DELETE;
import static com.rtg.position.output.Offset.INSERT;
import static com.rtg.position.output.Offset.TOTAL;
import static com.rtg.position.output.Offset.ZERO;

import com.rtg.util.StringUtils;
import com.rtg.util.format.FormatInteger;
import com.rtg.util.format.FormatReal;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 * Used to explore the distribution of gapped values
 * for long read Ngs matches.
 * Currently fixed for KAREN defaults of sub prob = 0.001 indel prob= 0.009
 * (these were the numbers from 454).
 */
public class GappedDistribution extends IntegralAbstract {

  /**
   * Check if p is a valid probability.
   * @param p value to be checked.
   * @return true (exception thrown is not valid).
   * @throws IllegalArgumentException when p is invalid.
   */
  public static boolean prob(final double p) {
    if (p <= 0.0 || p > 1.0 || Double.isNaN(p) || Double.isInfinite(p)) {
      throw new IllegalArgumentException(String.valueOf(p));
    }
    return true;
  }

  /*
   * Get a default version for use with testing etc.
   */
  static PositionDistributionParams distrParams(final int maxGap) {
    return new PositionDistributionParams(0.001, 0.009, maxGap, 0);
  }

  /** Probability of an equal match. */
  private final double mPequal;

  /** Probability of a substitution. */
  private final double mPsubs;

  /** Probability of an equal or a substitution. */
  private final double mPeqs;

  /** Total probability of either an insertion or a deletion. */
  private final double mPiord;

  /** Probability of an individual insertion or a deletion. */
  private final double mPindel;

  private final int mStepSize;

  private final int mWordSize;

  private final int mGap;

  private final int mMult;

  private final int mILength;

  private final int mJLength;

  private final double[][][] mDistribution;

  /** Probability that a whole word will match (no substitutions or indels). */
  private final double mMatch;

  /**
   * Prefered way to construct.
   * @param params parameters that specify everything.
   */
  public GappedDistribution(final PositionParams params) {
    this(params.build().stepSize(), params.build().stepSize(), params.output().distribution());
  }


  /**
   * For testing only.
   * @param stepSize step size.
   * @param wordSize word size.
   * @param params parameters that determine largest gap and probabilities.
   */
  GappedDistribution(final int stepSize, final int wordSize, final PositionDistributionParams params) {
    mStepSize = stepSize;
    mWordSize = wordSize;
    mGap = params.maxGap();
    mMult = (mGap + mStepSize - 1) / mStepSize + 1;
    mILength = stepSize * mMult + 1;
    mJLength = stepSize * (mMult + 1) + 1;
    mDistribution = new double[mILength][][];
    for (int i = 0; i < mDistribution.length; i++) {
      mDistribution[i] = new double[mJLength][];
      for (int j = 0; j < mDistribution[i].length; j++) {
        mDistribution[i][j] = new double[Offset.values().length];
      }
    }
    mPsubs = params.subsProb();
    mPiord = params.indelOpenProb();
    //System.err.println("subs=" + mPsubs + " indel=" + mPiord);
    mPequal = 1.0 - mPsubs - mPiord;
    mPindel = mPiord / 2.0;
    mPeqs = mPequal + mPsubs;
    mMatch = Math.pow(mPequal, mStepSize);
    distribute();
  }

  private void distribute() {
    mDistribution[0][0][ZERO.ordinal()] = 1.0;
    mDistribution[0][0][TOTAL.ordinal()] = 1.0;

    //traverse the array diagonally - k == i + j
    final int kLength = mILength + mJLength - 1;
    for (int k = 1; k < kLength; k++) {
      for (int i = Math.max(0, k - mJLength + 1); i < Math.min(mILength, k + 1); i++) {
        final boolean boundary = i > 0 && (i % mStepSize) == 0;
        final int j = k - i;
        //System.err.println("k=" + k + " i=" + i + " j=" + j + " il=" + mILength + " jl=" + mJLength);
        final double[] dist = mDistribution[i][j];
        //diagonal
        if (i >= 1 && j >= 1) {
          final double[] dd = mDistribution[i - 1][j - 1];
          dist[ZERO.ordinal()] = (dd[ZERO.ordinal()] + dd[INSERT.ordinal()] + dd[DELETE.ordinal()]) * mPeqs;
        }
        //horizontal
        if (j >= 1) {
          final double[] dd = mDistribution[i][j - 1];
          final double s = dd[ZERO.ordinal()] + dd[INSERT.ordinal()] * 2.0;
          dist[INSERT.ordinal()] = s * mPindel;
        }
        //vertical
        if (i >= 1) {
          final double[] dd = mDistribution[i - 1][j];
          dist[DELETE.ordinal()] = (dd[ZERO.ordinal()] + dd[DELETE.ordinal()] * 2.0) * mPindel;
        }
        //match case
        if (boundary && j >= mStepSize) {
          assert i >= mStepSize;
          final int is = i - mStepSize;
          final int js = j - mStepSize;
          final double[] dz = mDistribution[is][js];
          final double match = dz[ZERO.ordinal()] * mMatch;
          dist[ZERO.ordinal()] -= match;
        }
        //Total
        dist[TOTAL.ordinal()] = dist[ZERO.ordinal()] + dist[INSERT.ordinal()] + dist[DELETE.ordinal()];
      }
    }
  }

  /**
   * Compute the offset that the values will have when word size is not a multiple of the step size.
   * In a separate method so that it can be tested.
   * @param stepSize size of steps between hits.
   * @param wordSize size of word.
   * @return the offset (&gt;= 0, &lt; step size).
   */
  static int stepOffset(final int stepSize, final int wordSize) {
    assert stepSize <= wordSize;
    final int mod = wordSize % stepSize;
    if (mod == 0) {
      return 0;
    }
    final int x = stepSize - mod;
    assert x >= 0 && x < stepSize : x;
    return x;
  }

  /**
   * Construct a gap probabilities table.
   *
   * @return a valid gap probabilities table.
   */
  GapScorer probabilities() {
    final int jlen = mGap + 1;
    assert jlen <= mJLength;
    //System.err.println("ILength=" + mILength + " step=" + mStepSize + " isLen=" + isLen);
    //assert (isLen - 1) * mStep < mILength : "gap=" + mGap + " step=" + mStep + " isLen=" + isLen + " mILength=" + mILength;
    final double[][] probs = new double[mILength][];
    for (int i = 0; i < probs.length; i++) {
      probs[i] = new double[jlen];
    }
    //compute the threshold
    double max = 0.0;
    //System.err.println("j=[" + mGap + ".." + mJLength + ")");
    //System.err.println("i=[" + offset + ".." + mILength + "+" + mStepSize + ")");
    for (int j = mGap + 1; j < mJLength; j++) {
      for (int i = mGap + 1; i < mILength; i++) {
        final double d = mDistribution[i][j][TOTAL.ordinal()];
        //System.err.println("  max=" + max + " d=" + d + " i=" + i + " j=" + j);
        if (d > max) {
          //System.err.println("**max=" + max + " d=" + d + " i=" + i + " j=" + j);
          max = d;
        }
      }
    }
    final double min = mDistribution[mGap][mGap][TOTAL.ordinal()];
    final double gapThreshold = max > min ? min : max;
    for (int j = 0; j <= mGap; j++) {
      for (int i = 0, is = 0; is < mILength; i++, is++) {
        final double d = mDistribution[i][j][TOTAL.ordinal()];
        probs[is][j] = d < gapThreshold ? 0.0 : d;
      }
    }

    //compute max probabilities for use with READ scoring
    final double[] maxd = new double[mILength];
    for (int i = 0; i < mILength; i++) {
      double mx = 0.0;
      for (int j = 0; j < mJLength; j++) {
        final double d = mDistribution[i][j][TOTAL.ordinal()];
        if (d > mx) {
          mx = d;
        }
      }
      assert mx >= 0.0 && mx <= 1.0 && !Double.isNaN(mx);
      maxd[i] = mx;
    }
    //System.err.println(Arrays.toString(maxd));
    final GapProbabilitiesScorer gapProbabilities = new GapProbabilitiesScorer(gapThreshold, probs, maxd);
    gapProbabilities.globalIntegrity();
    return gapProbabilities;
  }

  private static final FormatInteger INDEX_FORMAT = new FormatInteger(3);
  private static final FormatReal PROB_FORMAT = new FormatReal(2, 4);

  @Override
  public void toString(final StringBuilder out) {
    for (int i = mILength - 1; i >= 0; i--) {
      for (int l = 0; l < Offset.values().length; l++) {
        switch (l) {
        case 0: out.append("[").append(INDEX_FORMAT.format(i)).append("]=");
        break;
        case 1: out.append("     i");
        break;
        case 2: out.append("     d");
        break;
        case 3: out.append("     t");
        break;
        default: throw new IllegalStateException();
        }
        double total = 0.0;
        for (int j = 0; j < mJLength; j++) {
          out.append(PROB_FORMAT.format(mDistribution[i][j][l])).append(i == j ? "=" : " ");
          total += mDistribution[i][j][l];
        }
        out.append("*").append(PROB_FORMAT.format(total));
        out.append(StringUtils.LS);
      }
      out.append(StringUtils.LS);
    }
  }

  @Override
  public boolean integrity() {
    Exam.assertTrue(prob(mPequal));
    Exam.assertTrue(prob(mPsubs));
    Exam.assertTrue(prob(mPeqs));
    Exam.assertTrue(prob(mPiord));
    Exam.assertTrue(prob(mPindel));
    Exam.assertEquals(1.0, mPequal + mPsubs + mPiord, 1e-8);
    Exam.assertTrue(mStepSize >= 1);
    Exam.assertTrue(mStepSize <= mWordSize);
    Exam.assertTrue(mMult >= 1);
    Exam.assertTrue(mILength > mStepSize);
    Exam.assertTrue(mJLength > mStepSize);
    Exam.assertEquals(mILength, mStepSize * mMult + 1);
    Exam.assertEquals(mJLength, mStepSize * (mMult + 1) + 1);
    return true;
  }

  /**
   * @param args command line arguments - args[0] step size &gt;= 1.
   */
  public static void main(final String[] args) {
    final int step = Integer.parseInt(args[0]);
    if (step < 1) {
      throw new IllegalArgumentException(step + "");
    }
    final int word = Integer.parseInt(args[1]);
    if (word < 1) {
      throw new IllegalArgumentException(word + "");
    }
    final int mult = Integer.parseInt(args[2]);
    if (mult < 1) {
      throw new IllegalArgumentException(mult + "");
    }
    final GappedDistribution gd = new GappedDistribution(step, word, distrParams(mult));
    System.out.print(gd.toString());
  }
}
