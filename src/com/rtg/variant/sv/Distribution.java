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

package com.rtg.variant.sv;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.util.StringUtils;
import com.rtg.util.Utils;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 */
@TestClass(value = {"com.rtg.variant.sv.DistributionArrayTest"})
public abstract class Distribution extends IntegralAbstract {

  protected final int mLo;

  protected final int mHi;

  /**
   * @param lo low index of the distribution (inclusive).
   * @param hi high index of the distribution (exclusive).
   */
  public Distribution(int lo, int hi) {
    super();
    mLo = lo;
    mHi = hi;
    //System.err.println("lo=" + lo + " hi=" + hi);
  }

  /**
   * Get the value in the distribution.
   * @param index position in the distribution (from -radius to radius-1 inclusive).
   * @return the value.
   */
  public double get(final int index) {
    if (index < mLo || index >= mHi) {
      throw new IndexOutOfBoundsException("index=" + index + " lo=" + mLo + " hi=" + mHi);
    }
    return getValue(index);
  }

  /**
   * Return a log-based signal for this distribution with the given sam counts
   * @param counts the sam counts
   * @param label a textual label
   * @return the signal
   */
  public Signal getSignalLn(SamCounts counts, String label) {
    return new SignalDistributionLn(counts, this, label);
  }

  /**
   * Get a value in the distribution.
   * This version is local and does not check the index values.
   * @param index position in the distribution (0 based).
   * @return the value.
   */
  protected abstract double getValue(int index);

  /**
   * Lowest index in the distribution (inclusive).
   * @return the low index.
   */
  int lo() {
    return mLo;
  }

  /**
   * Highest index in the distribution (exclusive).
   * @return the high index.
   */
  int hi() {
    return mHi;
  }

  String dump() {
    final StringBuilder sb = new StringBuilder();
    for (int i = lo(); i < hi(); ++i) {
      sb.append(i).append(" ").append(Utils.realFormat(get(i), 4)).append(StringUtils.LS);
    }
    return sb.toString();
  }

  @Override
  public boolean integrity() {
    //System.err.println(mLo + ":" + mHi);
    Exam.assertTrue("" + mLo, mLo < 0);
    Exam.assertTrue(mHi > 0);
    return true;
  }

  @Override
  public boolean globalIntegrity() {
    integrity();
    for (int i = lo(); i < hi(); ++i) {
      final double v = getValue(i);
      Exam.assertTrue("" + v, !Double.isNaN(v) && !Double.isInfinite(v));
    }
    return true;
  }

}
