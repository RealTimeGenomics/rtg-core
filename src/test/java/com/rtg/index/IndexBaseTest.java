/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
    for (int i = 0; i < x.length; ++i) {
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
