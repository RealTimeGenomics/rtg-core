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

import static com.rtg.reader.CgUtils.CG2_EXPECTED_OVERLAP;
import static com.rtg.reader.CgUtils.CG2_OVERLAP_POSITION;
import static com.rtg.reader.CgUtils.CG2_RAW_READ_LENGTH;
import static com.rtg.reader.CgUtils.CG_EXPECTED_LARGE_GAP;
import static com.rtg.reader.CgUtils.CG_EXPECTED_OVERLAP;
import static com.rtg.reader.CgUtils.CG_OVERLAP_POSITION;
import static com.rtg.reader.CgUtils.CG_RAW_READ_LENGTH;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.util.diagnostic.SpyCounter;
import com.rtg.util.integrity.Exam;
import com.rtg.util.machine.MachineType;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;


/**
 * Calculates a banded matrix of Complete Genomics alignment probabilities.
 *
 */
@TestClass({"com.rtg.variant.realign.ScoreMatrixCGTestSuite"})
public class ScoreMatrixCG extends ScoreMatrix {

  static final int OVERLAP_GAP = RealignParamsImplementation.CG_OVERLAP;
  static final int SMALL_GAP = RealignParamsImplementation.CG_SMALL_GAP;
  static final int LARGE_GAP = RealignParamsImplementation.CG_LARGE_GAP;
  static final int OVERLAP2_GAP = RealignParamsImplementation.CG_OVERLAP2;

  private static final SpyCounter SPY = new SpyCounter("ScoreMatrixCG");
  private final double[][] mGapDistributionsPoss;
  private final boolean mVersion1;

  /**
   * @param arith helper object that does the arithmetic so that this code can be independent of the representation.
   * @param params CG priors and gap probabilities
   */
  public ScoreMatrixCG(final PossibilityArithmetic arith, final RealignParams params) {
    super(arith, params);
    mGapDistributionsPoss = params.gapDistributionPoss(arith);
    mVersion1 = params.machineType() == MachineType.COMPLETE_GENOMICS;
    SPY.increment();
  }

  private double gapFreqPoss(final int gap, final int width) {
    return mGapDistributionsPoss[gap][width - mParams.gapStart(gap)];
  }

  @Override
  public void setEnv(Environment env) {
    super.setEnv(env);
    if (mLength != (mVersion1 ? CG_RAW_READ_LENGTH : CG2_RAW_READ_LENGTH)) {
      throw new IllegalArgumentException("CG read is not the right length was: " + mLength);
    }
  }

  protected static int rowOffsetCGV1(final int row, final int maxShift) {
    final int diagonal = row - maxShift - 1;
    if (row <= CG_OVERLAP_POSITION) {
      return diagonal;
    } else if (row <= 25) {
      return diagonal - CG_EXPECTED_OVERLAP;
    } else {
      return diagonal - CG_EXPECTED_OVERLAP + CG_EXPECTED_LARGE_GAP;
    }
  }

  protected static int rowOffsetCGV2(final int row, final int maxShift) {
    final int diagonal = row - maxShift - 1;
    if (row <= CG2_OVERLAP_POSITION) {
      return diagonal;
    } else {
      return diagonal - CG2_EXPECTED_OVERLAP;
    }
  }

  @Override
  protected int rowOffset(final int row) {
    return mVersion1 ? rowOffsetCGV1(row, mEnv.maxShift()) : rowOffsetCGV2(row, mEnv.maxShift());
  }

  @Override
  protected boolean isFragmentStart(final int row) {
    if (mVersion1) {
      return row == 6 || row == 16 || row == 26;
    } else {
      return row == 11 || row == 20;
    }
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

    for (int j = 0; j <= lastCol; ++j) {
      // delete: we sum the probabilities over all gap sizes
      double del = mArith.zero();
      for (int gapSize = gapStart; gapSize <= gapEnd; ++gapSize) {
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
      for (int gapSize = gapStart; gapSize <= gapEnd; ++gapSize) {
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
   * The CG version does more complex calculations after each gap.
   */
  @Override
  protected void calculateProbabilities() {
    if (mVersion1) {
      calculateProbabilitiesV1();
    } else {
      calculateProbabilitiesV2();
    }
  }

  private void calculateProbabilitiesV1() {
    int row = 0;
    calculateInitialRow(row++, mDeleteStartPoss, mMatchStartPoss);
    while (row <= CG_OVERLAP_POSITION) {
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

  private void calculateProbabilitiesV2() {
    int row = 0;
    calculateInitialRow(row++, mDeleteStartPoss, mMatchStartPoss);
    while (row <= CG2_OVERLAP_POSITION) {
      calculateRow(row++);
    }
    calculateCGGap(row++, OVERLAP2_GAP);
    while (row <= mLength) {
      calculateRow(row++);
    }
    calculateEnd();
  }

  @Override
  public boolean integrity() {
    final boolean result = super.integrity();
    for (int gap = 0; gap < 2; ++gap) {
      // check that we don't have any columns with no numbers in them
      Exam.assertTrue(Math.abs(mParams.gapStart(gap)) < mEnv.maxShift());
      Exam.assertTrue(Math.abs(mParams.gapEnd(gap)) < mEnv.maxShift());
    }
    return result;
  }

}
