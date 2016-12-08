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
public class CodeDiploidTest extends TestCase {

  public void testHaploid() {
    final CodeDiploid tc = new CodeDiploid(4);
    assertEquals(4, tc.rangeSize());
    for (int i = 0; i < 4; ++i) {
      assertTrue(tc.homozygous(i));
      assertEquals(tc.a(i), tc.bc(i));
    }
    for (int i = 4; i < 10; ++i) {
      assertFalse(tc.homozygous(i));
      assertFalse(tc.a(i) == tc.bc(i));
    }

  }

  public void testValid() {
    final CodeDiploid tc = new CodeDiploid(4);
    assertFalse(tc.valid(-1));
    assertFalse(tc.valid(10));
    for (int i = 0; i < 10; ++i) {
      assertTrue(tc.valid(i));
    }
  }

  private void check(final CodeDiploid tc, final int n, final int i, final int j) {
    assertEquals(j, tc.j(n));
    assertEquals(i, tc.i(n));
    assertEquals(i + j, tc.k(n));
  }

  public void test4() {
    final CodeDiploid tc = new CodeDiploid(4);
    assertEquals(10, tc.size());
    check(tc, 0, 0, 0);
    check(tc, 1, 1, 0);
    check(tc, 2, 0, 1);
    check(tc, 3, 2, 0);
    check(tc, 4, 1, 1);
    check(tc, 5, 0, 2);
    check(tc, 6, 3, 0);
    check(tc, 7, 2, 1);
    check(tc, 8, 1, 2);
    check(tc, 9, 0, 3);
  }

  public void test2() {
    final CodeDiploid tc = new CodeDiploid(2);
    assertEquals(3, tc.size());
    check(tc, 0, 0, 0);
    check(tc, 1, 1, 0);
    check(tc, 2, 0, 1);
  }

  private void checkab(final Code tc, final int n, final int a, final int b) {
    assertTrue(b <= a);
    assertEquals(a, tc.a(n));
    assertEquals(b, tc.b(n));
  }

  public void test4ab() {
    final CodeDiploid tc = new CodeDiploid(4);
    assertEquals(10, tc.size());
    checkab(tc, 0, 0, 0);
    checkab(tc, 1, 1, 1);
    checkab(tc, 2, 2, 2);
    checkab(tc, 3, 3, 3);
    checkab(tc, 4, 1, 0);
    checkab(tc, 5, 2, 1);
    checkab(tc, 6, 3, 2);
    checkab(tc, 7, 2, 0);
    checkab(tc, 8, 3, 1);
    checkab(tc, 9, 3, 0);
  }

  public void testCode1() {
    final CodeDiploid tc = new CodeDiploid(4);
    for (int i = 0; i < 4; ++i) {
      assertEquals(i, tc.code(i));
      assertEquals(i, tc.code(i, i));
    }
  }

  public void testCode2() {
    final CodeDiploid tc = new CodeDiploid(4);
    for (int i = 0; i < tc.size(); ++i) {
      final int a = tc.a(i);
      final int b = tc.homozygous(i) ? a : tc.b(i);
      assertEquals(i, tc.code(a, b));
      assertEquals(i, tc.code(b, a));
    }
  }
}
