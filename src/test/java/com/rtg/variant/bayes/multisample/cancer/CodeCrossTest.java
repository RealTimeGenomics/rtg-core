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

package com.rtg.variant.bayes.multisample.cancer;

import com.rtg.variant.bayes.Code;

import junit.framework.TestCase;

/**
 */
public class CodeCrossTest extends TestCase {

  public void test() {

  }

  public void testValid() {
    final Code tc = new CodeCross(4);
    assertEquals(16, tc.size());
    assertFalse(tc.valid(-1));
    assertFalse(tc.valid(16));
    for (int i = 0; i < 16; ++i) {
      assertTrue(tc.valid(i));
    }
  }

  private void checkab(final Code tc, final int n, final int a, final int b) {
    assertEquals(a, tc.a(n));
    assertEquals(b, tc.b(n));
    assertEquals(b, tc.bc(n));
  }

  public void test4ab() {
    final Code tc = new CodeCross(3);
    assertEquals(9, tc.size());
    checkab(tc, 0, 0, 0);
    checkab(tc, 1, 0, 1);
    checkab(tc, 2, 0, 2);
    checkab(tc, 3, 1, 0);
    checkab(tc, 4, 1, 1);
    checkab(tc, 5, 1, 2);
    checkab(tc, 6, 2, 0);
    checkab(tc, 7, 2, 1);
    checkab(tc, 8, 2, 2);
  }

  public void testCode2() {
    final Code tc = new CodeCross(4);
    assertEquals(4, tc.rangeSize());
    for (int i = 0; i < tc.size(); ++i) {
      final int a = tc.a(i);
      assertTrue(0 <= a && a < 4);
      final int b = tc.bc(i);
      assertTrue(0 <= b && b < 4);
      assertEquals(i, tc.code(a, b));
    }
  }

  public void testCode() {
    final Code tc = new CodeCross(3);
    for (int i = 0; i < 3; ++i) {
      assertEquals(i, tc.code(i));
      assertTrue(tc.homozygous(tc.code(i, i)));
      for (int j = 0; j < 3; ++j) {
        if (i == j) {
          continue;
        }
        assertFalse(tc.homozygous(tc.code(i, j)));
      }

    }
  }
}
