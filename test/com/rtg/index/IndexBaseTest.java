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

import com.rtg.util.TestUtils;
import com.rtg.util.array.ArrayUtils;
import com.rtg.util.array.CommonIndex;

/**
 */
public class IndexBaseTest extends IndexSimpleTest {

  private CommonIndex dynamic(final int[] x) {
    final CommonIndex a = com.rtg.util.array.intindex.IntCreate.createIndex(x.length);
    for (int i = 0; i < x.length; i++) {
      a.set(i, x[i]);
    }
    return a;
  }

  public void testIndexState() {
    TestUtils.testEnum(IndexBase.IndexState.class, "[PRE_ADD, ADD, FROZEN]");
  }
  public void testIsSortedLong() {
    assertTrue(ArrayUtils.isSorted(dynamic(new int[0]), 0));
    assertTrue(ArrayUtils.isSorted(dynamic(new int[]{5}), 0));
    assertTrue(ArrayUtils.isSorted(dynamic(new int[]{5}), 1));
    assertTrue(ArrayUtils.isSorted(dynamic(new int[]{5, 6}), 0));
    assertTrue(ArrayUtils.isSorted(dynamic(new int[]{5, 6}), 1));
    assertTrue(ArrayUtils.isSorted(dynamic(new int[]{5, 6}), 2));
    assertTrue(ArrayUtils.isSorted(dynamic(new int[]{5, 3}), 0));
    assertTrue(ArrayUtils.isSorted(dynamic(new int[]{5, 3}), 1));
    assertFalse(ArrayUtils.isSorted(dynamic(new int[]{5, 3}), 2));
    assertTrue(ArrayUtils.isSorted(dynamic(new int[]{5, 5}), 2));
  }
  public void testIsSortedSegment() {
    assertTrue(ArrayUtils.isSorted(dynamic(new int[0]), 0, 0));
    assertTrue(ArrayUtils.isSorted(dynamic(new int[]{5, 4, 3}), 0, 1));
    assertTrue(ArrayUtils.isSorted(dynamic(new int[]{5, 4, 3}), 2, 3));
    assertTrue(ArrayUtils.isSorted(dynamic(new int[]{5, 4, 3}), 2, 3));
    assertTrue(ArrayUtils.isSorted(dynamic(new int[]{3, 5, 4, 5, 3, 5}), 0, 2));
    assertTrue(ArrayUtils.isSorted(dynamic(new int[]{3, 5, 4, 5, 3, 5}), 2, 4));
    assertTrue(ArrayUtils.isSorted(dynamic(new int[]{3, 5, 4, 5, 3, 4, 5}), 4, 7));
    assertFalse(ArrayUtils.isSorted(dynamic(new int[]{3, 5, 4, 5, 3, 5}), 0, 3));
    assertFalse(ArrayUtils.isSorted(dynamic(new int[]{3, 5, 4, 5, 3, 5}), 1, 4));
  }
  public void testIsSortedStrictLong() {
    assertTrue(ArrayUtils.isSortedStrict(dynamic(new int[0]), 0));
    assertTrue(ArrayUtils.isSortedStrict(dynamic(new int[]{5}), 0));
    assertTrue(ArrayUtils.isSortedStrict(dynamic(new int[]{5}), 1));
    assertTrue(ArrayUtils.isSortedStrict(dynamic(new int[]{5, 6}), 0));
    assertTrue(ArrayUtils.isSortedStrict(dynamic(new int[]{5, 6}), 1));
    assertTrue(ArrayUtils.isSortedStrict(dynamic(new int[]{5, 6}), 2));
    assertTrue(ArrayUtils.isSortedStrict(dynamic(new int[]{5, 3}), 0));
    assertTrue(ArrayUtils.isSortedStrict(dynamic(new int[]{5, 3}), 1));
    assertFalse(ArrayUtils.isSortedStrict(dynamic(new int[]{5, 3}), 2));
    assertFalse(ArrayUtils.isSortedStrict(dynamic(new int[]{5, 5}), 2));
  }


}
