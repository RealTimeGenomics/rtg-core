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
package com.rtg.util.array.zeroindex;

import static com.rtg.util.StringUtils.LS;

import junit.framework.TestCase;

/**
 */
public class ZeroIndexTest extends TestCase {

  public void testLengthNegative() {
    checkError(-1, "length=-1");
  }

  public void testLengthEmpty() {
    final ZeroIndex index = ZeroCreate.createIndex(0);
    // test toString() when all zeroes
    assertEquals("Index [0]" + LS, index.toString());
    assertEquals(0, index.bytes());
    assertEquals("     ", index.formatValue().blanks());
  }

  public void test() {
    final ZeroIndex index = ZeroCreate.createIndex(2);
    index.set(1, 0);
    index.setSigned(1, 0);
    try {
      index.get(2);
      fail();
    } catch (final IndexOutOfBoundsException e) {
      assertEquals("Index out of bounds:2 : 2", e.getMessage());
    }
    try {
      index.get(-1);
      fail();
    } catch (final IndexOutOfBoundsException e) {
      assertEquals("Index out of bounds:-1 : 2", e.getMessage());
    }
    assertEquals(0, index.get(0));
    assertEquals(0, index.getSigned(0));
    assertEquals("Index [2]" + LS, index.toString());
  }

  public void testConstant() {
    final ZeroIndex index = new ZeroIndex(2, -42);
    assertEquals(-42, index.get(0));
    assertEquals(-42, index.getSigned(0));
    assertEquals("Index [2]" + LS + "-42 constant" + LS, index.toString());
  }

  public void checkError(final long length, final String expected) {
    try {
      ZeroCreate.createIndex(length);
      fail("expected exception: " + expected);
    } catch (final Exception e) {
      assertEquals(expected, e.getMessage());
    }
  }

  public void testExtend() {
    final ZeroIndex index = ZeroCreate.createIndex(2);
    assertEquals(2, index.length());
    assertEquals(2, index.extendBy(3));
    assertEquals(5, index.length());
  }

  public void testTrim() {
    final ZeroIndex index = ZeroCreate.createIndex(5);
    assertEquals(5, index.length());
    index.trim(3);
    assertEquals(3, index.length());
  }
}
