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
package com.rtg.variant.cnv;

import java.io.IOException;
import java.util.Arrays;

import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 */
public class LeastSquaresModel extends IntegralAbstract {

  private static final int INITIAL_SIZE = 1000;
  private static final double SPLIT_CONSTANT = 0.5 * (Math.log(2 * Math.PI) + 1.0);

  private final String mSeqName;
  private final LsmOutput[] mOut;
  private final double mExtraPenalty;
  private final int mMinBlocks;
  //private final double mCorrection;

  private double[] mSum;
  private double[] mSum2;

  private int mNext = 1;
  private final boolean mExtraPenaltyOff;

  /**
   * Constructor
   * @param seqName unique id for sequence
   * @param minBlocks minimum number of blocks in each side of a split
   * @param extraPenalty extra amount a split needs to overcome to be considered an improvement
   * @param extraPenaltyOff flag if extra penalty is switched off
   * @param correction the bucket size
   * @param out where to write output
   */
  public LeastSquaresModel(final String seqName, final int minBlocks, final double correction, final double extraPenalty,
      final boolean extraPenaltyOff, final LsmOutput... out) {
    mSeqName = seqName;
    mExtraPenalty = extraPenalty;
    mExtraPenaltyOff = extraPenaltyOff;
    mMinBlocks = minBlocks;
    //mCorrection = 1 + Math.ceil(correction);
    mOut = out;
    mSum = new double[INITIAL_SIZE];
    mSum2 = new double[INITIAL_SIZE];
  }

  private void rescale() {
    assert mNext == mSum.length && mNext == mSum2.length;
    mSum = rescale(mSum, mNext);
    mSum2 = rescale(mSum2, mNext);
  }

  private double[] rescale(final double[] array, final int length) {
    return Arrays.copyOf(array, (length + 1) * 2);
  }

  /**
   * Add a value to the LeastSquaresModel
   * @param value the number to add
   */
  public void add(final double value) {
    if (mNext == mSum.length) {
      rescale();
    }
    mSum[mNext] = mSum[mNext - 1] + value;
    mSum2[mNext] = mSum2[mNext - 1] + value * value;
    ++mNext;
    //System.err.println("add next=" + mNext);
  }

  double mGlobalVariance = 0;

  /**
   * Scan a region from start to end inclusive (1 based which the arrays are not).
   * @throws IOException when writing output
   */
  public void fit() throws IOException {
    mGlobalVariance = variance(1, mNext);
    scan(1, mNext);
  }

  private double variance(final int start, final int end) {
    final int n = end - start;
    assert n > 1;
    final double d1 = mSum[end - 1] - mSum[start - 1];
    final double x1 = d1 * d1 / n;
    final double x2 = mSum2[end - 1] - mSum2[start - 1];
    return (x2 - x1) / (n - 1);
  }

  private double score(final int start, final int end) {
    final int n = end - start;
    final double v = variance(start, end);
    if (v <= 0) {
      return n * SPLIT_CONSTANT;
    }
    //    final double s = -0.5 * n * Math.log(v);
    return -0.5 * n * Math.log(1 + v);
  }

  private double splitPenalty(final int start, final int end) {
    // Represents the cost of encoding the position of the break point, should be ~log(n)
    final int n = end - start;
    //    return 10 * Math.log(n);
    //    return Math.log(n) + 10000 / n;

    // global variance always 0??
    return Math.log(n); // + (mCorrection / n) * mGlobalVariance;
    //    return Math.log(n); // John's original
  }

  /**
   * Scan a region from start to end inclusive (1 based which the arrays are not).
   * @param start position starting at 1.
   * @param end position (exclusive)
   * @throws IOException when writing output.
   */
  public void scan(final int start, final int end) throws IOException {
    if (start < 1 || start > end) {
      throw new IllegalArgumentException("start=" + start + " end=" + end);
    }
    //best is highest scoring split.
    int bestPosition = 0;
    if (end - start > mMinBlocks) {
      double best = mExtraPenaltyOff ? score(start, end) + splitPenalty(start, end)
          : score(start, end) + splitPenalty(start, end) + mExtraPenalty;
      //System.err.println("scan start=" + start + " end=" + end + " best=" + best);

      for (int i = start + mMinBlocks; i < end - mMinBlocks; ++i) {
        final double x1 = score(start, i);
        final double x2 = score(i, end);
        final double x = x1 + x2;
        //System.err.println("x=" + x + " i=" + i);
        if (x > best) {
          best = x;
          bestPosition = i;
        }
      }
    }
    if (bestPosition == 0) {
      //no split is better write this region
      for (final LsmOutput element : mOut) {
        final double sum = mSum[end - 1] - mSum[start - 1];
        final double sum2 = mSum2[end - 1] - mSum2[start - 1];
        element.out(mSeqName, start, end, sum, sum2);
      }
    } else {
      scan(start, bestPosition);
      scan(bestPosition, end);
    }
  }

  /**
   * Close any files or other resources.
   * @throws IOException IO Exceptions occur sometimes
   */
  public void close() throws IOException {
    for (final LsmOutput out : mOut) {
      out.close();
    }
  }

  @Override
  public void toString(final StringBuilder sb) {
    sb.append("LeastSquaresModel");
  }

  @Override
  public boolean globalIntegrity() {
    integrity();
    for (int i = 1; i < mNext; ++i) {
      Exam.assertTrue(mSum2[i] >= 0.0);
      Exam.assertTrue(mSum2[i] >= mSum2[i - 1]);
    }
    return true;
  }

  @Override
  public boolean integrity() {
    Exam.assertTrue(mNext >= 1);
    Exam.assertTrue(mMinBlocks > 0);
    Exam.assertEquals(mSum.length, mSum2.length);
    Exam.assertTrue(mSum.length >= mNext);
    Exam.assertEquals(0.0, mSum[0]);
    Exam.assertEquals(0.0, mSum2[0]);
    return true;
  }


}
