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

import com.rtg.util.diagnostic.SpyCounter;
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

  private static final SpyCounter SPY = new SpyCounter("ScoreMatrix");

  /**
   * A score matrix with the given maximum band width.
   * @param arith helper object that does the arithmetic so that this code can be independent of the representation.
   * @param params the machine error model and related parameters.
   */
  protected ScoreMatrix(PossibilityArithmetic arith, final RealignParams params) {
    super(arith, params);
    SPY.increment();
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
    for (int j = 0; j < mWidth; j++) {
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
    for (int j = mWidth - 1; j >= 0; j--) {
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
    for (int i = 1; i <= mLength; i++) {
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
