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

import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 * Store a one dimensional interval. Can be oriented either up or down.
 * a is inclusive, b is exclusive.
 */
final class Interval extends IntegralAbstract {
  private final int mX;
  private final int mY;

  Interval(int a, int b) {
    mX = a;
    mY = b;
  }

  Interval negative() {
    return new Interval(-mX, -mY);
  }

  /**
   * Get a.
   * @return Returns the a.
   */
  public int getA() {
    return mX;
  }

  /**
   * Get b.
   * @return Returns the b.
   */
  public int getB() {
    return mY;
  }

  @Override
  public boolean equals(Object obj) {
    if (obj == null) {
      return false;
    }
    if (obj == this) {
      return true;
    }
    final Interval that = (Interval) obj;

    return this.getA() == that.getA() && this.getB() == that.getB();
  }

  @Override
  public int hashCode() {
    // pacify findbugs
    return super.hashCode();
  }

  @Override
  public boolean integrity() {
    Exam.assertFalse(mX == mY);
    return true;
  }
}
