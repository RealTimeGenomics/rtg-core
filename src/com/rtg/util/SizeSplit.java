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
package com.rtg.util;

/**
 * Facilities for splitting a range as evenly as possible.
 * Handy for doing things like allocating work to threads.
 */
public class SizeSplit {
  private final int mN;

  private final int mD;

  private final int mS;

  private final int mR;

  /**
   * @param n number of unit to be divided.
   * @param d the number of groups they are to be divided into.
   */
  public SizeSplit(final int n, final int d) {
    assert n >= 0 && d >= 1;
    mN = n;
    mD = d;
    mS = n / d;
    mR = n - d * mS;
    assert mR >= 0;
    assert mS >= 0;
  }

  /**
   * Get the first position (zero based, inclusive) allocated to "unit" i.
   * @param i unit of work ranging from 0 to d inclusive.
   * @return the last allocated value (zero based, inclusive).
   */
  public int start(final int i) {
    assert 0 <= i && i <= mD;
    return i <= mR ? i * (mS + 1) : mR + i * mS;
  }

  @Override
  public String toString() {
    return "SizeSplit: n=" + mN + " d=" + mD + " s=" + mS + " r=" + mR;
  }

}
