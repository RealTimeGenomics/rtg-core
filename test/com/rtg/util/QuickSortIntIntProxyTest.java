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

import java.util.Arrays;

import junit.framework.TestCase;

/**
 */
public class QuickSortIntIntProxyTest extends TestCase {
  public void test() {
    final int[] moreints = {50, 42, 10, 45, 20, 20};
    final int[] ints = {0, 1, 2, 3, 4, 5};
    final QuickSortIntIntProxy proxy = new QuickSortIntIntProxy(moreints, ints, false);
    QuickSort.sort(proxy);
    final int[] expected = {0, 3, 1, 4, 5, 2};
    assertTrue("Expected: " + Arrays.toString(expected) + " actual: " + Arrays.toString(ints), Arrays.equals(expected, ints));
    assertTrue(Arrays.equals(new int[]{50, 45, 42, 20, 20, 10}, moreints));

    final QuickSortIntIntProxy proxy2 = new QuickSortIntIntProxy(moreints, ints, true);
    QuickSort.sort(proxy2);
    assertTrue(Arrays.equals(new int[]{10, 20, 20, 42, 45, 50}, moreints));
  }
}
