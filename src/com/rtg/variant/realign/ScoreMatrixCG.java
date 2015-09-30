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

import static com.rtg.reader.CgUtils.CG_RAW_READ_LENGTH;

import com.rtg.util.diagnostic.SpyCounter;
import com.rtg.util.integrity.Exam;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;


/**
 * Calculates a banded matrix of Complete Genomics alignment probabilities.
 *
 */
public class ScoreMatrixCG extends ScoreMatrix {

  static final int OVERLAP_GAP = RealignParamsImplementation.CG_OVERLAP;
  static final int SMALL_GAP = RealignParamsImplementation.CG_SMALL_GAP;
  static final int LARGE_GAP = RealignParamsImplementation.CG_LARGE_GAP;

  private static final SpyCounter SPY = new SpyCounter("ScoreMatrixCG");
  private final double[][] mGapDistributionsPoss;

  /**
   * @param arith helper object that does the arithmetic so that this code can be independent of the representation.
   * @param params CG priors and gap probabilities
   */
  public ScoreMatrixCG(final PossibilityArithmetic arith, final RealignParams params) {
    super(arith, params);
    mGapDistributionsPoss = params.gapDistributionPoss(arith);
    SPY.increment();
  }

  protected final double gapFreqPoss(final int gap, final int width) {
    return mGapDistributionsPoss[gap][width - mParams.gapStart(gap)];
  }

  @Override
  public void setEnv(Environment env) {
    super.setEnv(env);
    if (mLength != CG_RAW_READ_LENGTH) {
      throw new IllegalArgumentException("CG read is not the right length was: " + mLength);
    }
  }

  protected static int rowOffsetCG(final int row, final int maxShift) {
    final int diagonal = row - maxShift - 1;
    if (row <= 5) {
      return diagonal;
    } else if (row <= 25) {
      return diagonal - 2;
    } else {
      return diagonal - 2 + 6;
    }
  }

  @Override
  protected int rowOffset(final int row) {
    return rowOffsetCG(row, mEnv.maxShift());
  }

  /**
   * This is similar to <code>calculateDelete</code>, but always uses the
   * open-delete penalty, never the extend-delete penalty.
   * It should be used whenever a CG gap is crossed, even if the gap is zero,
   * because the chemistry is discontinuous at each gap.
   * @param i row
   * @param j column
   * @return natural log.
   */
  protected double calculateOpenDelete(final int i, final int j) {
    final double s = mArith.add(delete(i, j), match(i, j));
    //System.err.println("mDeleteOpenInFourPoss=" + mDeleteOpenInFourPoss + " s=" + s + " delete=" + delete(i, j) + " match=" + match(i, j));
    return mArith.multiply(s, mDeleteOpenInFourPoss);
  }

  protected double calculateMatchCGGap(final int i, final int j) {
    final double del = mArith.multiply(delete(i, j), mOneMinusDeleteOpenPoss);
    final double match = mArith.multiply(match(i, j), mOneMinusDeleteInsertOpenPoss);
    final double ins = mArith.multiply(insert(i, j), mOneMinusInsertExtendPoss);
    return mArith.add(mArith.add(del, match), ins);
  }

  /**
   * Calculates a CG row after a variable size gap.
   * @param row the first row after the gap.
   * @param whichGap which gap to handle (see <code>RealignParams</code>).
   */
  protected void calculateCGGap(final int row, final int whichGap) {
    final int gapStart = mParams.gapStart(whichGap);
    final int gapEnd = mParams.gapEnd(whichGap);
    //    final byte re = mEnv.read(row - 1);
    final int lastCol = mWidth - 1;
    final int thisRowOffset = rowOffset(row);
    final int offset = thisRowOffset - rowOffset(row - 1);

    for (int j = 0; j <= lastCol; j++) {
      // delete: we sum the probabilities over all gap sizes
      double del = mArith.zero();
      for (int gapSize = gapStart; gapSize <= gapEnd; gapSize++) {
        final int prevCol = j - (gapSize - offset);
        if (0 <= prevCol && prevCol < mWidth) {
          final double from = calculateOpenDelete(row - 1, prevCol);
          final double gapFreq = gapFreqPoss(whichGap, gapSize);
          del = mArith.add(del, mArith.multiply(from, gapFreq));
        }
      }
      setDelete(row, j, del);

      // match or mismatch: we sum the probabilities over all gap sizes
      double mm = mArith.zero();
      for (int gapSize = gapStart; gapSize <= gapEnd; gapSize++) {
        final int prevCol = j - (gapSize - offset) - 1;
        //        final byte te = mEnv.template(thisRowOffset + j);
        if (0 <= prevCol && prevCol < mWidth) {
          final double gapFreq = gapFreqPoss(whichGap, gapSize);
          final double s1 = mArith.multiply(calculateMatchCGGap(row - 1, prevCol), gapFreq);
          mm = mArith.add(mm, s1);
        }
      }
      setMatch(row, j, mArith.multiply(mm, matchEq(row - 1, j)));

      // insert is the same as usual
      setInsert(row, j, j == 0 ? mArith.zero() : calculateInsert(row, j - 1));
    }
  }

  /**
   * The CG version does more complex calculations after each gap (5, 15, 25).
   */
  @Override
  protected void calculateProbabilities() {
    int row = 0;
    calculateInitialRow(row++, mDeleteStartPoss, mMatchStartPoss);
    while (row <= 5) {
      calculateRow(row++);
    }
    calculateCGGap(row++, OVERLAP_GAP);
    while (row <= 15) {
      calculateRow(row++);
    }
    calculateCGGap(row++, SMALL_GAP);
    while (row <= 25) {
      calculateRow(row++);
    }
    calculateCGGap(row++, LARGE_GAP);
    while (row <= mLength) {
      calculateRow(row++);
    }
    calculateEnd();
  }

  @Override
  public boolean integrity() {
    final boolean result = super.integrity();
    for (int gap = 0; gap < 2; gap++) {
      // check that we don't have any columns with no numbers in them
      Exam.assertTrue(Math.abs(mParams.gapStart(gap)) < mEnv.maxShift());
      Exam.assertTrue(Math.abs(mParams.gapEnd(gap)) < mEnv.maxShift());
    }
    return result;
  }

}
