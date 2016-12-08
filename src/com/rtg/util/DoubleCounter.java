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

import com.reeltwo.jumble.annotations.TestClass;

/**
 * A mutable counter.
 */
@TestClass("com.rtg.util.DoubleMultiSetTest")
public class DoubleCounter {

  private double mCount = 0;

  /**
   * Get the current count
   * @return the current count
   */
  public double count() {
    return mCount;
  }

  /**
   * Increment the counter.
   */
  public void increment() {
    ++mCount;
  }

  /**
   * Increment the counter.
   * @param count the amount to increment by.
   */
  public void increment(double count) {
    assert count >= 0;
    mCount += count;
  }

}
