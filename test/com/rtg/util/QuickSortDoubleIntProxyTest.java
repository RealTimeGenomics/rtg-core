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
public class QuickSortDoubleIntProxyTest extends TestCase {
  public void test() {
    final double[] doubles = {5.0, 4.2, 1.0, 4.5, 2.0, 2.0};
    final int[] ints = {0, 1, 2, 3, 4, 5};
    final QuickSortDoubleIntProxy proxy = new QuickSortDoubleIntProxy(doubles, ints);
    QuickSort.sort(proxy);
    final int[] expected = {0, 3, 1, 4, 5, 2};
    assertTrue("Expected: " + Arrays.toString(expected) + " actual: " + Arrays.toString(ints), Arrays.equals(expected, ints));
    assertTrue(Arrays.equals(new double[]{5.0, 4.5, 4.2, 2.0, 2.0, 1.0}, doubles));
  }
}
