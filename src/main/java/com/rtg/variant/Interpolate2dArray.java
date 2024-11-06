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
