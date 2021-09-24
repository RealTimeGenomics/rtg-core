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

package com.rtg.variant.sv.discord;


/**
 */
class BreakpointPosition {

  private final int mLo;
  private final int mPosition;
  private final int mPositionY;
  private final int mHi;

  BreakpointPosition(int lo, int position, int hi, int positiony) {
    if (lo > position || hi < position) {
      throw new IllegalArgumentException();
    }
    mLo = lo;
    mPosition = position;
    mPositionY = positiony;
    mHi = hi;
  }

  int lo() {
    return mLo;
  }

  int position() {
    return mPosition;
  }

  int positionAlt() {
    return mPositionY;
  }

  int hi() {
    return mHi;
  }

  @Override
  public String toString() {
    return "BreakpointPosition [" + "lo=" + mLo + ", position=" + mPosition + ", hi=" + mHi + ", positionY=" + mPositionY + "]";
  }
}
