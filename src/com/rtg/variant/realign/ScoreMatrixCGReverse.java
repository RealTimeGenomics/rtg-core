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

import com.rtg.variant.util.arithmetic.PossibilityArithmetic;


/**
 * Calculates a banded matrix of Complete Genomics alignment probabilities
 * in the reverse direction (from the end of the read, back to the beginning).
 *
 */
public class ScoreMatrixCGReverse extends ScoreMatrixReverse {

  static final int OVERLAP_GAP = 0;
  static final int SMALL_GAP = 1;
  static final int LARGE_GAP = 2;

  private final double[][] mGapDistributionsPoss;

  /**
   * @param arith helper object that does the arithmetic so that this code can be independent of the representation.
   * @param params CG priors and gap probabilities
   */
  public ScoreMatrixCGReverse(final PossibilityArithmetic arith, final RealignParams params) {
    super(arith, params);
    mGapDistributionsPoss = params.gapDistributionPoss(arith);
  }

  private double gapFreqPoss(final int gap, final int width) {
    return mGapDistributionsPoss[gap][width - mParams.gapStart(gap)];
  }

  @Override
  protected int rowOffset(final int row) {
    return ScoreMatrixCG.rowOffsetCGV1(row, mEnv.maxShift());
  }

  /**
   * Calculates a CG row before a variable size gap.
   * @param row the row just before the gap.
   * @param whichGap which gap to handle (see <code>RealignParams</code>).
   */
  protected void calculateCGGap(final int row, final int whichGap) {
    assert 0 < row && row < mLength;
    final int gapStart = mParams.gapStart(whichGap);
    final int gapEnd = mParams.gapEnd(whichGap);
    final int offset = rowOffset(row + 1) - rowOffset(row);

    for (int j = mWidth - 1; j >= 0; j--) {
      double dd = mArith.zero();
      double mm = mArith.zero();
      double ii = mArith.zero();
      for (int gapSize = gapStart; gapSize <= gapEnd; gapSize++) {
        final int matchCol = j - offset + gapSize + 1;
        if (0 <= matchCol && matchCol < mWidth) {
          final double gapLn = gapFreqPoss(whichGap, gapSize);
          final double ma = matchEq(row, matchCol);

          // delete
          final double m1 = mArith.multiply(gapLn, ma);
          final double m2 = mArith.multiply(m1, mOneMinusDeleteOpenPoss);
          final double m3 = mArith.multiply(m2, match(row + 1, matchCol));
          dd = mArith.add(dd, m3);
          final double m4 = mArith.multiply(gapLn, mDeleteOpenPoss);
          if (matchCol > 0) {
            final double m5 = mArith.multiply(delete(row + 1, matchCol - 1), mOneInFourPoss);
            final double m = mArith.multiply(m4, m5);
            dd = mArith.add(dd, m);
          }

          // match
          final double m6 = mArith.multiply(m1, mOneMinusDeleteInsertOpenPoss);
          final double m7 = mArith.multiply(m6, match(row + 1, matchCol));
          mm = mArith.add(mm, m7);
          if (matchCol > 0) {
            final double m8 = mArith.multiply(m4, mOneInFourPoss);
            final double m9 = mArith.multiply(m8, delete(row + 1, matchCol - 1));
            mm = mArith.add(mm, m9);
          }

          // insert
          final double mm1 = mArith.multiply(m1, mOneMinusInsertExtendPoss);
          final double mm2 = mArith.multiply(mm1, match(row + 1, matchCol));
          ii = mArith.add(ii, mm2);
        }
      }
      // now add in the (horizontal) insert links, which are independent of the gap size.
      if (j < mWidth - 1) {
        final double mm1 = mArith.multiply(mInsertOpenPoss, insert(row, j + 1));
        mm = mArith.add(mm, mm1);
        final double mm2 = mArith.multiply(mInsertExtendPoss, insert(row, j + 1));
        ii = mArith.add(ii, mm2);
      }
      setDelete(row, j,  dd);
      setMatch(row, j,  mm);
      setInsert(row, j,  ii);
    }
  }

  /**
   * The CG version does more complex calculations after each gap (5, 15, 25).
   */
  @Override
  protected void calculateProbabilities() {
    int row = mLength;
    final double one = mArith.one();
    calculateInitialRow(row--, one, one);
    while (row > 25) {
      calculateRow(row--);
    }
    calculateCGGap(row--, LARGE_GAP);
    while (row > 15) {
      calculateRow(row--);
    }
    calculateCGGap(row--, SMALL_GAP);
    while (row > 5) {
      calculateRow(row--);
    }
    calculateCGGap(row--, OVERLAP_GAP);
    while (row >= 0) {
      calculateRow(row--);
    }
    calculateEnd();
  }
}
