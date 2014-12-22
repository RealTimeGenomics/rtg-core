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

package com.rtg.variant.bayes;

import junit.framework.TestCase;

/**
 */
public class CodeHaploidTest extends TestCase {

  public void testHaploid() {
    final Code tc = new CodeHaploid(4);
    assertEquals(4, tc.rangeSize());
    assertEquals(4, tc.size());
    for (int i = 0; i < 4; i++) {
      assertTrue(tc.homozygous(i));
      assertEquals(tc.a(i), tc.bc(i));
    }
  }

  public void testValid() {
    final Code tc = new CodeHaploid(4);
    assertFalse(tc.valid(-1));
    assertFalse(tc.valid(4));
    for (int i = 0; i < 4; i++) {
      assertTrue(tc.valid(i));
      assertEquals(i, tc.a(i));
      assertEquals(i, tc.code(i));
      assertEquals(i, tc.code(i, i));
    }

    try {
      tc.b(0);
      fail();
    } catch (final UnsupportedOperationException e) {
      //expected
    }
    try {
      tc.code(-1);
      fail();
    } catch (final IllegalArgumentException e) {
      assertEquals("-1", e.getMessage());
    }
    try {
      tc.code(1, 2);
      fail();
    } catch (final UnsupportedOperationException e) {
      //expected
    }
  }
}
