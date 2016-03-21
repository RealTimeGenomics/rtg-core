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

package com.rtg.variant;

/**
 * Performs interpolation within a 2d array
 */
class Interpolate2dArray implements Interpolate {
  private final int mMinX;
  private final int mMinY;
  private final int mStepX;
  private final int mStepY;
  private final int[][] mCurve;

  Interpolate2dArray(int[][] curve, int minX, int minY, int stepX, int stepY) {
    mMinX = minX;
    mMinY = minY;
    mStepX = stepX;
    mStepY = stepY;
    mCurve = curve;
  }

  @Override
  public int getValue(int pos) {
    return mCurve[mMinX + pos * mStepX][mMinY + pos * mStepY];
  }

  @Override
  public void setValue(int pos, int value) {
    mCurve[mMinX + pos * mStepX][mMinY + pos * mStepY] = value;
  }

  @Override
  public int minPos() {
    return 0;
  }

  @Override
  public boolean inBounds(int pos) {
    final int x = mMinX + pos * mStepX;
    final int y = mMinY + pos * mStepY;
    return x >= 0 && x < mCurve.length && y >= 0 && y < mCurve[0].length;
  }

  @Override
  public boolean isMissing(int pos) {
    return getValue(pos) == -1;
  }

  /**
   * Return an interpolation down a column
   * @param curve curve to interpolate
   * @param column column of the curve to interpolate
   * @return a new Interpolate2dArray configured for the {@code column} column of {@code curve}
   */
  public static Interpolate2dArray column(int[][] curve, int column) {
    return new Interpolate2dArray(curve, 0, column, 1, 0);
  }
  /**
   * Return an interpolation across a row
   * @param curve curve to interpolate
   * @param row column of the curve to interpolate
   * @return a new Interpolate2dArray configured for the {@code row} row of {@code curve}
   */
  public static Interpolate2dArray row(int[][] curve, int row) {
    return new Interpolate2dArray(curve, row, 0, 0, 1);
  }
}
