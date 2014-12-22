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
package com.rtg.util.array;


import junit.framework.TestCase;

/**
 */
public class WrappedIntArrayTest extends TestCase {

  /**
   * Test method for {@link WrappedIntArray}.
   */
  public final void test() {
    final int[] a = {1, 2, 3};
    final WrappedIntArray wa = new WrappedIntArray(a);
    wa.integrity();
    assertEquals(3, wa.length());
    assertEquals("[1, 2, 3]", wa.toString());
    assertEquals(1, wa.get(0));
    assertEquals(2, wa.get(1));
    assertEquals(3, wa.get(2));
    try {
      wa.get(-1);
      fail();
    } catch (final ArrayIndexOutOfBoundsException e) {
      //expected
    }
    try {
      wa.get(3);
      fail();
    } catch (final ArrayIndexOutOfBoundsException e) {
      //expected
    }
  }

  public final void testLong() {
    final long[] a = {1, 2, 3};
    final WrappedIntArray wa = new WrappedIntArray(a);
    wa.integrity();
    assertEquals(3, wa.length());
    assertEquals("[1, 2, 3]", wa.toString());
    assertEquals(1, wa.get(0));
    assertEquals(2, wa.get(1));
    assertEquals(3, wa.get(2));
  }

  public final void testPaired() {
    final int[] a = {1, 2, 3};
    final int[] b = {4, 5, 6};
    final WrappedIntArray wa = new WrappedIntArray(a, b);
    wa.integrity();
    assertEquals(6, wa.length());
    assertEquals("[1, 4, 2, 5, 3, 6]", wa.toString());
    assertEquals(1, wa.get(0));
    assertEquals(4, wa.get(1));
    assertEquals(2, wa.get(2));
    assertEquals(5, wa.get(3));
    assertEquals(3, wa.get(4));
    assertEquals(6, wa.get(5));
  }
}
