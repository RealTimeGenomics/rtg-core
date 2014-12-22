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
package com.rtg.util.array.objectindex;


import junit.framework.TestCase;

/**
 * Tests for <code>ObjectIndex</code>.
 *
 */
public abstract class AbstractObjectIndexTest extends TestCase {

  protected static final int STEP = 1000;
  protected static final int LENGTH = 1000007; //this isnt a multiple of two because then it might not catch many tricky cases

  /** Local new line convention */
  private static final String LS = System.lineSeparator();

  /**
   * Constructor (needed for JUnit)
   * @param name A string which names the tests.
   */
  public AbstractObjectIndexTest(final String name) {
    super(name);
  }

  protected abstract ObjectIndex<Integer> create(final long length);

  protected abstract ObjectIndex<Integer> create(final long length, final int bits);

  public void testMasks() {
    assertEquals(ObjectIndex.INT_MASK & ObjectIndex.HIGH_MASK, 0);
    assertEquals(ObjectIndex.HIGH_SIGNED_MASK & ObjectIndex.HIGH_MASK, ObjectIndex.HIGH_MASK);
    assertEquals(ObjectIndex.HIGH_SIGNED_MASK & ObjectIndex.INT_MASK, 1L << 31);
  }

  private static final String TO_STR = ""
    + "Index [100]" + LS
    + "[         0] " + null + ", 1, 2, " + null + ", " + null + ", " + null + ", " + null + ", " + null + ", " + null + ", " + null + "" + LS
    + "[        10] " + null + ", " + null + ", 12, " + null + ", " + null + ", " + null + ", " + null + ", " + null + ", " + null + ", " + null + "" + LS
    + "[        50] " + null + ", " + null + ", 52, " + null + ", " + null + ", " + null + ", " + null + ", " + null + ", " + null + ", " + null + "" + LS;

  public void testToString() {
    final ObjectIndex<Integer> index = create(100);
    index.set(1, 1);
    index.set(2, 2);
    index.set(12, 12);
    index.set(52, 52);
    final String str = index.toString();
    assertEquals(TO_STR, str);
  }

  private static final String TO_STR15 = ""
    + "Index [15]" + LS
    + "[         0] " + null + ", 1, 2, " + null + ", " + null + ", " + null + ", " + null + ", " + null + ", " + null + ", " + null + "" + LS
    + "[        10] " + null + ", " + null + ", 12, " + null + ", " + null + "" + LS;

  public void testToString15() {
    final ObjectIndex<Integer> index = create(15);
    index.set(1, 1);
    index.set(2, 2);
    index.set(12, 12);
    final String str = index.toString();
    assertEquals(TO_STR15, str);
  }

  public void testShortToString() {
    final ObjectIndex<Integer> index = create(3);
    index.set(1, 1);
    final String str = index.toString();
    assertEquals("Index [3]" + LS + "[         0] " + null + ", 1, " + null + "" + LS, str);
  }

  public void testLength() {
    final int le = 42000;
    final ObjectIndex<Integer> a = create(le);
    a.integrity();
    assertEquals(le, a.length());
    assertEquals(4 * le, a.bytes());
    a.close();
    a.integrity();
  }

  public void testLength0() {
    final int le = 0;
    final ObjectIndex<Integer> a = create(le);
    a.integrity();
    assertEquals(le, a.length());
    assertEquals(4 * le, a.bytes());
    a.close();
    a.integrity();
  }

  public void testBadLength1() {
    try {
      create(Integer.MIN_VALUE);
      fail("Expected NegativeArraySizeException");
    } catch (final NegativeArraySizeException e) {
      // expected
    }
    try {
      create(-1);
      fail("Expected NegativeArraySizeException");
    } catch (final NegativeArraySizeException e) {
      // expected
    }
  }

  public void testIntensiveSetGet() {
    //needed for subtle errors in underlying mapping in disk backed cases
    final int length = 100; // > 2^5 (see ShortDiskChunksTest and ShortChunksTest - also not a multiple of 2^3
    final ObjectIndex<Integer> a = create(length, 3);
    a.integrity();
    assertEquals(length, a.length());
    for (int i = 0; i < a.length(); i++) {
      assertEquals(null, a.get(i));
      final int j = i * 3;
      a.set(i, j);
    }
    for (int i = 0; i < a.length(); i++) {
      final int j = i * 3;
      assertEquals(Integer.valueOf(j), a.get(i));
    }
    a.close();
    a.integrity();
  }

  public void testGetSetLong() {
    final ObjectIndex<Integer> a = create(LENGTH);
    a.integrity();
    assertEquals(LENGTH, a.length());
    for (int i = 0; i < a.length(); i += STEP) {
      assertEquals(null, a.get(i));
      final int j = i * 3;
      a.set(i, j);
    }
    for (int i = 0; i < a.length(); i += STEP) {
      final int j = i * 3;
      assertEquals(Integer.valueOf(j), a.get(i));
    }
    a.close();
    a.integrity();
    try {
      a.set(0L, 0);
      fail("Expected Exception");
    } catch (final RuntimeException e) {
      // expected
    }
    try {
      a.get(0L);
      fail("Expected Exception");
    } catch (final RuntimeException e) {
      // expected
    }
  }

  public void testSwap() {
    final ObjectIndex<Integer> a = create(LENGTH);
    a.integrity();
    assertEquals(LENGTH, a.length());
    for (int i = 0; i < a.length(); i += STEP) {
      assertEquals(null, a.get(i));
      final int j = i * 3;
      a.set(i, j);
      assertEquals(Integer.valueOf(j), a.get(i));
    }
    for (int i = 0; i < a.length() - 1; i += STEP) {
      final int j = i * 3;
      assertEquals(i + ":" + 0, null, a.get(i + 1));
      assertEquals(i + ":" + j, Integer.valueOf(j), a.get(i));
      a.swap(i, i + 1);
      assertEquals(i + ":" + j, Integer.valueOf(j), a.get(i + 1));
      assertEquals(i + ":" + 0, null, a.get(i));
      a.swap(i, i + 1);
      assertEquals(i + ":" + 0, null, a.get(i + 1));
      assertEquals(i + ":" + j, Integer.valueOf(j), a.get(i));
    }
    a.close();
    a.integrity();
    try {
      a.swap(0L, 0L);
      fail("Expected Exception");
    } catch (final RuntimeException e) {
      // expected
    }
  }


  public void testEdges() {
    final ObjectIndex<Integer> a = create(LENGTH);
    a.integrity();
    assertEquals(LENGTH, a.length());
    a.set(LENGTH - 1, 1);
    assertEquals(Integer.valueOf(1), a.get(LENGTH - 1));
    checkLimits(a, LENGTH);
    a.close();
    a.integrity();
    checkClose(a, LENGTH);
  }

  private void checkLimits(final ObjectIndex<Integer> a, final int length) {
    a.set(length - 1, 42);
    assertEquals(Integer.valueOf(42), a.get(length - 1));
    try {
      a.get(length);
      fail("Exception expected");
    } catch (final RuntimeException | AssertionError e) {
      //expected
    }
    try {
      a.set(length, 0);
      fail("Exception expected");
    } catch (final RuntimeException e) {
      //expected
    }
    a.set(0, 1);
    assertEquals(Integer.valueOf(1), a.get(0));
    try {
      a.get(-1);
      fail("Exception expected");
    } catch (final RuntimeException | AssertionError e) {
      //expected
    }
    try {
      a.set(-1, 0);
      fail("Exception expected");
    } catch (final RuntimeException | AssertionError e) {
      //expected
    }
  }

  private void checkClose(final ObjectIndex<Integer> a, final int length) {
    assertEquals(length, a.length());
    assertEquals(4L * length, a.bytes());
    try {
      a.set(length - 1, 42);
      fail("Exception expected");
    } catch (final RuntimeException e) {
      //expected
    }
    try {
      a.get(length - 1);
      fail("Exception expected");
    } catch (final RuntimeException e) {
      //expected
    }
    try {
      a.set(0, 1);
      fail("Exception expected");
    } catch (final RuntimeException e) {
      //expected
    }
    try {
      a.get(0);
      fail("Exception expected");
    } catch (final RuntimeException e) {
      //expected
    }
  }
}


