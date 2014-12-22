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

package com.rtg.sam;

import java.util.Set;
import java.util.TreeSet;

import junit.framework.TestCase;

/**
 */
public class RandomArrayListTest extends TestCase {

  public void test() {
    check(1, 2, 3, 4);
    check();
    check(1);
  }

  private void check(final Integer... values) {
    final RandomArrayList<Integer> ra = new RandomArrayList<>(42, values);
    final Set<Integer> set = new TreeSet<>();
    final Set<Integer> set0 = new TreeSet<>();
    for (Integer value : values) {
      final Integer next = ra.next();
      set.add(next);
      set0.add(value);
    }
    assertNull(ra.next());
    assertEquals(set0, set);
  }
}
