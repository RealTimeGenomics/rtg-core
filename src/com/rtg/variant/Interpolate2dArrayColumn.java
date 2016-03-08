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
 * @author kurt
 */
class Interpolate2dArrayColumn implements Interpolate {
  private final int mColumn;
  private final int[][] mCurve;

  public Interpolate2dArrayColumn(int column, int[][] curve) {
    mColumn = column;
    mCurve = curve;
  }

  @Override
  public int getValue(int pos) {
    return mCurve[pos][mColumn];
  }

  @Override
  public void setValue(int pos, int value) {
    mCurve[pos][mColumn] = value;
  }

  @Override
  public int minPos() {
    return 0;
  }

  @Override
  public int maxPos() {
    return mCurve.length;
  }
}
