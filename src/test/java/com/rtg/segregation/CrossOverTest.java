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

package com.rtg.segregation;

import junit.framework.TestCase;

/**
 */
public class CrossOverTest extends TestCase {

  public void test0() {
    final PatternArray pa = new PatternArray("0", "1");
    final CrossOver co = new CrossOver(3, false, 42, 2, 1, pa);
    co.integrity();
    assertTrue(pa == co.pattern());
    assertEquals("mo 42 2 1", co.toString());

    final int[] bedSearch0 = co.searchLevel(3, 2);
    assertEquals(3, bedSearch0[0]);
    assertEquals(1, bedSearch0[1]);

    final int[] bedSearch1 = co.searchLevel(3, 1);
    assertEquals(3, bedSearch1[0]);
    assertEquals(2, bedSearch1[1]);
  }

  public void test1() {
    final PatternArray pa = new PatternArray("0", "1");
    final CrossOver co = new CrossOver(3, true, 0, 1, 2, pa);
    co.integrity();
    assertTrue(pa == co.pattern());
    assertEquals("fa 0 1 2", co.toString());

    final int[] bedSearch0 = co.searchLevel(1, 3);
    assertEquals(2, bedSearch0[0]);
    assertEquals(3, bedSearch0[1]);

    final int[] bedSearch1 = co.searchLevel(2, 3);
    assertEquals(1, bedSearch1[0]);
    assertEquals(3, bedSearch1[1]);
  }
}
