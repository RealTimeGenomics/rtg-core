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
package com.rtg.index;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.util.array.CommonIndex;

/**
 * Function for swapping two array elements
 */
@TestClass(value = {"com.rtg.util.QuickSortTest"})
public class Swapper {

  /**
   * Swapper for two index arrays.
   * @param primary primary index array.
   * @param secondary secondary index array.
   */
  public Swapper(final CommonIndex primary, final CommonIndex secondary) {
    mPrimary = primary;
    mSecondary = secondary;
  }

  private final CommonIndex mPrimary;
    private final CommonIndex mSecondary;

  /**
   * Swap the two elements a and b.
   * @param a first element.
   * @param b second element.
   */
  public void swap(final long a, final long b) {
//    System.err.println("CI swap: " + a + " " + b);
    //System.err.println(a + ":" + b);
    final long x = mPrimary.get(a);
    mPrimary.set(a, mPrimary.get(b));
    mPrimary.set(b, x);
    final long y = mSecondary.get(a);
    mSecondary.set(a, mSecondary.get(b));
    mSecondary.set(b, y);
  }
}
