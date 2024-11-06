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

package com.rtg.variant.realign;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;


/**
 * Calculates a banded matrix of Complete Genomics alignment probabilities
 * in the reverse direction (from the end of the read, back to the beginning).
 *
 */
@TestClass({"com.rtg.variant.realign.ScoreMatrixCGReverseTestSuite"})
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

    for (int j = mWidth - 1; j >= 0; --j) {
      double dd = mArith.zero();
      double mm = mArith.zero();
      double ii = mArith.zero();
      for (int gapSize = gapStart; gapSize <= gapEnd; ++gapSize) {
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
