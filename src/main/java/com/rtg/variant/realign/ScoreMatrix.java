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

import com.rtg.variant.util.arithmetic.PossibilityArithmetic;


/**
 * Calculates a banded matrix of alignment probabilities.
 * This is similar to the <code>Gotoh</code> algorithm, but calculates the
 * probability of all paths to a given point, rather than just the best path.
 *
 * This does just the forward direction matrix.
 *
 */
public class ScoreMatrix extends AbstractAllPaths {

  /**
   * A score matrix with the given maximum band width.
   * @param arith helper object that does the arithmetic so that this code can be independent of the representation.
   * @param params the machine error model and related parameters.
   */
  protected ScoreMatrix(PossibilityArithmetic arith, final RealignParams params) {
    super(arith, params);
  }

  protected final double calculateDelete(final int i, final int j) {
    if (j == mWidth) {
      return mZeroPoss;
    }
    final double sum = mArith.add(
        mArith.multiply(mDeleteExtendPoss, delete(i, j)),
        mArith.multiply(mDeleteOpenPoss, match(i, j)));
    return mArith.multiply(sum, mOneInFourPoss);
  }

  /**
   * Calculate the forward diagonal probabilities, starting from cell <code>(i,j)</code>.
   * This does NOT include the match/mismatch probability of the destination cell.
   *
   * @param i zero-based read position
   * @param j column
   * @return Possibility
   */
  protected final double calculateMatch(final int i, final int j) {
    final double del =   mArith.multiply(delete(i, j), mOneMinusDeleteExtendPoss);
    final double match = mArith.multiply(match(i, j), mOneMinusDeleteInsertOpenPoss);
    final double ins =   mArith.multiply(insert(i, j), mOneMinusInsertExtendPoss);
    return mArith.add(mArith.add(del, match), ins);
  }

  protected final double calculateInsert(final int i, final int j) {
    if (i == mLength || j < 0) {
      return mZeroPoss;
    }
    return mArith.add(
        mArith.multiply(mInsertExtendPoss, insert(i, j)),
        mArith.multiply(mInsertOpenPoss, match(i, j)));
  }

  protected final void calculateRow(final int i) {
    for (int j = 0; j < mWidth; ++j) {
      setDelete(i, j, calculateDelete(i - 1, j + 1));
      matchIt(i, j);
      setInsert(i, j, calculateInsert(i, j - 1));
    }
  }

  protected void matchIt(final int i, int j) {
    setMatch(i, j, mArith.multiply(calculateMatch(i - 1, j), matchEq(i - 1, j)));
  }

  /**
   * Calculates <code>mEndScores</code>
   */
  protected final void calculateEnd() {
    double sumLn = mZeroPoss;
    for (int j = mWidth - 1; j >= 0; --j) {
      // last insert row is all zero probabilities for most matrices, but is not initialised in the delta matrices
      // assert mInsert[endRow][j] == mZeroPoss;
      sumLn = mArith.add(sumLn, mArith.add(delete(mLength, j), match(mLength, j)));
      mEndScores[j] = sumLn;
    }
  }

  @Override
  protected void calculateProbabilities() {
    calculateInitialRow(0, mDeleteStartPoss, mMatchStartPoss);
    //insert in read and delete on template
    for (int i = 1; i <= mLength; ++i) {
      calculateRow(i);
    }
    calculateEnd();
  }

  @Override
  public final double totalScoreLn() {
    return mArith.poss2Ln(mEndScores[0]);
  }

  @Override
  public final double totalScore() {
    return mArith.poss2Prob(mEndScores[0]);
  }

  @Override
  public double total() {
    return mEndScores[0];
  }

  @Override
  public final boolean underflow() {
    return mArith.underflow(mEndScores[0]);
  }

  /**
   * Calculates the likelihood that the end of the read is at or after
   * the given template position.
   *
   * @param templatePos an absolute template position.
   * @return the natural log of the probability.
   */
  public final double readEndsAfterLn(final int templatePos) {
    final int col = templatePos - mEnv.absoluteTemplatePosition(0) - rowOffset(mLength);
    final double resultPoss;
    if (col < 0) {
      resultPoss = mOnePoss;  // the end is definitely to the right of templatePos
    } else if (col >= mWidth) {
      resultPoss = mZeroPoss;  // end is definitely to the left of templatePos
    } else {
      resultPoss = mArith.divide(mEndScores[col], mEndScores[0]);
    }
    return mArith.poss2Ln(resultPoss);
  }
}
