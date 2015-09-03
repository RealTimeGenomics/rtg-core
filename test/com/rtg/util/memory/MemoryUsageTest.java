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
package com.rtg.util.memory;

import java.lang.ref.SoftReference;

import junit.framework.TestCase;

/**
 */
public class MemoryUsageTest extends TestCase {

  private MemoryUsage mMem = null;


  static class Foo {
    protected int mD1, mD2;
  }

  public void testInteger() {
    tstObject(0, MemoryUsage.INT_SIZE);
  }

  public void testFoo() {
    tstObject(new MemoryUsageFoo(), 2 * MemoryUsage.INT_SIZE);
  }

  void tstObject(final Object obj, final int isize) {
    final int size = MemoryUsage.roundUp(MemoryUsage.OBJECT_SIZE + isize);
    mMem.calculateUsage(obj);
    assertEquals("mem stats\n" + mMem, size, mMem.getSize());
    assertEquals("mem stats\n" + mMem, size, mMem.getHardSize());
    assertEquals("mem stats\n" + mMem, 0, mMem.getSoftSize());
    mMem.calculateUsage(obj, new Object[]{obj});
    assertEquals("mem stats\n" + mMem, 0, mMem.getSize());
    assertEquals("mem stats\n" + mMem, 0, mMem.getHardSize());
    assertEquals("mem stats\n" + mMem, 0, mMem.getSoftSize());

  }

  public void testIntegerArray() {
    final Integer[] ia = new Integer[]
      {0, 1, 2, 3};
    final int asize = MemoryUsage.roundUp(MemoryUsage.ARRAY_SIZE
      + 4 * MemoryUsage.REF_SIZE);
    final int isize = 4 * MemoryUsage.roundUp(MemoryUsage.OBJECT_SIZE
      + MemoryUsage.INT_SIZE);
    tst(ia, asize + isize);
  }

  public void testintArray() {
    final int[] ia = new int[4];
    final int asize = MemoryUsage.roundUp(MemoryUsage.ARRAY_SIZE
      + 4 * MemoryUsage.INT_ARRAY);
    tst(ia, asize);
  }


  /** Test when same object appears more than once.  */
  public void testIntegerDuplicates() {
    final Integer[] ia = new Integer[4];
    ia[0] = 0;
    ia[1] = ia[0];
    ia[2] = 2;
    ia[3] = ia[2];

    final int asize = MemoryUsage.roundUp(MemoryUsage.ARRAY_SIZE
      + 4 * MemoryUsage.REF_SIZE);
    final int isize = 2 * MemoryUsage.roundUp(MemoryUsage.OBJECT_SIZE
      + MemoryUsage.INT_SIZE);
    tst(ia, asize + isize);
  }

  public void testIllegalState() {
    try {
      mMem.getHardIterator();
      fail("Expected IllegalStateException");
    } catch (final IllegalStateException ex) {
      // Expected
    }
    try {
      mMem.getSoftIterator();
      fail("Expected IllegalStateException");
    } catch (final IllegalStateException ex) {
      // Expected
    }
    try {
      mMem.getHardSize();
      fail("Expected IllegalStateException");
    } catch (final IllegalStateException ex) {
      // Expected
    }
    try {
      mMem.getSoftSize();
      fail("Expected IllegalStateException");
    } catch (final IllegalStateException ex) {
      // Expected
    }
    try {
      mMem.getSize();
      fail("Expected IllegalStateException");
    } catch (final IllegalStateException ex) {
      // Expected
    }

  }


  /**
   * Test soft references.
   */
  public void testSoft() {
    final Object[] ia = new Object[6];
    final Integer i1 = 1;
    ia[0] = new SoftReference<>(i1);
    ia[1] = i1;
    final Integer i2 = 2;
    ia[2] = new SoftReference<>(i2);
    ia[3] = new SoftReference<>(i2);
    final SoftReference<?> s = new SoftReference<>(3);
    ia[4] = s;
    ia[5] = s;

    //array size - 6 references
    final int asize = MemoryUsage.roundUp(MemoryUsage.ARRAY_SIZE
      + 6 * MemoryUsage.REF_SIZE);
    //3 distinct integer objects
    final int is = MemoryUsage.roundUp(MemoryUsage.OBJECT_SIZE
      + MemoryUsage.INT_SIZE);
    final int isize = 3 * is;
    //4 distinct SoftReference objects
    final ClassMemory.Info info = ClassMemory.getMemoryInfo(SoftReference.class);
    final int rs = info.getSize();
    final int rsize = 4 * rs;
    final int size = asize + isize + rsize;
    //System.err.println("starting refs test");
    mMem.calculateUsage(ia);
    //System.err.println(mMem);
    assertTrue("Total size incorrect: " + mMem, mMem.getSize() >= size);
    assertEquals("mem stats\n" + mMem, asize + is, mMem.getHardSize());
    assertTrue("Soft size incorrect:" + mMem, mMem.getSoftSize() >= 4 * rs + 2 * is);
  }


  void tst(final Object obj, final int size) {
    mMem.calculateUsage(obj);
    assertEquals("mem stats\n" + mMem, size, mMem.getSize());
    assertEquals("mem stats\n" + mMem, size, mMem.getHardSize());
    assertEquals("mem stats\n" + mMem, 0, mMem.getSoftSize());
    mMem.calculateUsage(obj, new Object[]{obj});
    assertEquals("mem stats\n" + mMem, 0, mMem.getSize());
    assertEquals("mem stats\n" + mMem, 0, mMem.getHardSize());
    assertEquals("mem stats\n" + mMem, 0, mMem.getSoftSize());

  }

  /** Used by JUnit (called before each test method)  */
  @Override
  public void setUp() {
    mMem = new MemoryUsage();
  }


  /** Used by JUnit (called after each test method)  */
  @Override
  public void tearDown() {
    mMem = null;
  }

  public void testTotalAndClasses() {
    final Object o = new Object();
    final MemoryUsage.TotalAndClasses t = new MemoryUsageImpl(o).getTotalAndClasses();
    assertEquals(0, t.memSize(null));
    final String s = "hi";
    assertTrue(t.memSize(s) > 0);
    assertEquals(0, t.memSize(s));
    assertTrue(t.memSize(new Object()) > 0);
    assertEquals(0, t.memSize(o));
    assertTrue(t.memSize((float) 42) > 0);
    assertEquals(0, t.memSize(42));
  }

  public void testCount() {
    final MemoryUsage.Count c = new MemoryUsage.Count();
    assertEquals(0, c.getBytes());
    assertEquals(0, c.getReturned());
    assertEquals(0, c.getCount());
    c.push();
    c.pop();
    try {
      c.pop();
      fail();
    } catch (final IllegalStateException e) {
      assertEquals("Too many pops in memory calculation.", e.getMessage());
    }
    c.increment(42, 43);
    assertEquals(1, c.getCount());
    assertEquals(42, c.getBytes());
    assertEquals(43, c.getReturned());
  }

  public void testToString() {
    final MemoryUsage mu = new MemoryUsage();
    assertEquals("No object analysed.", mu.toString().trim());
    try {
      mu.getSize();
      fail();
    } catch (final IllegalStateException e) {
      assertEquals("No object has been analysed", e.getMessage());
    }
    mu.calculateUsage("x", new Object[0]);
    final String s = mu.toString();
    assertTrue(s.contains("MemoryUsage for object:"));
    assertTrue(s.contains("Total size:"));
    assertTrue(s.contains("Hard size:"));
    assertTrue(s.contains("Soft size:"));
    assertTrue(s.contains("Hard class statistics"));
    assertTrue(s.contains("Soft class statistics"));
    assertTrue(s.contains("Count"));
    assertTrue(s.contains("Size"));
    assertTrue(s.contains("Average"));
    assertTrue(s.contains("100.00%"));
    assertTrue(s.contains(" 0.00%"));
    assertTrue(s.contains("B  "));
    assertTrue(s.contains("#"));
    assertTrue(s.contains("Cumulative"));
    assertTrue(s.contains("Average"));
    assertTrue(s.contains("End MemoryUsage"));
    assertTrue(mu.getSize() > 0);
  }

  public void testConstantSizes() {
    assertEquals(4, MemoryUsage.refSize(short.class));
    assertEquals(4, MemoryUsage.refSize(char.class));
    assertEquals(4, MemoryUsage.refSize(int.class));
    assertEquals(4, MemoryUsage.refSize(float.class));
    assertEquals(8, MemoryUsage.refSize(long.class));
    assertEquals(8, MemoryUsage.refSize(double.class));
    assertEquals(4, MemoryUsage.refSize(byte.class));
    assertEquals(4, MemoryUsage.refSize(boolean.class));
    assertEquals(8, MemoryUsage.refSize(Object.class));
  }

  private class MemoryUsageImpl extends MemoryUsage {

    private final Object mO;

    public MemoryUsageImpl(final Object o) {
      this.mO = o;
      calculateUsage(new Object(), new Object[]{mO}, new Class<?>[]{Integer.class});
      mHardrefs = new IdentitySet();
    }

    TotalAndClasses getTotalAndClasses() {
      return new TotalAndClasses(false);
    }
  }
}

