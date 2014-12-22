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

import junit.framework.TestCase;

/**
 */
public class SizeSplitTest extends TestCase {

  public void test() {
    final SizeSplit ss = new SizeSplit(16, 3);
    assertEquals("SizeSplit: n=16 d=3 s=5 r=1", ss.toString());
    assertEquals(0, ss.start(0));
    assertEquals(6, ss.start(1));
    assertEquals(11, ss.start(2));
    assertEquals(16, ss.start(3));
  }

  public void test1() {
    check(16, 3);
    check(0, 100);
    check(100, 100);
    check(100, 99);
    check(100, 101);
  }

  void check(final int n, final int d) {
    final SizeSplit ss = new SizeSplit(n, d);
    final int s = n / d;
    int last = 0;
    for (int i = 1; i <= d; i++) {
      final int st = ss.start(i);
      final int diff = st - ss.start(i - 1);
      assertTrue(last <= st);
      assertTrue(diff == s || diff == s + 1);
      last = st;
    }
    assertEquals(0, ss.start(0));
    assertEquals(n, ss.start(d));
  }
}
