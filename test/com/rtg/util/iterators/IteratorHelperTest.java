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

package com.rtg.util.iterators;

import java.util.Iterator;

import junit.framework.TestCase;

/**
 */
public class IteratorHelperTest extends TestCase {

  private static class MockIterator extends IteratorHelper<Integer> {
    private final int mThreshold;
    private final int[] mX;
    private int mIndex = 0;

    MockIterator(int threshold, int... x) {
      mThreshold = threshold;
      mX = x;
    }

    @Override
    protected void step() {
      mIndex++;
    }

    @Override
    protected boolean isOK() {
      return mX[mIndex] >= mThreshold;
    }

    @Override
    protected boolean atEnd() {
      return mIndex >= mX.length;
    }

    @Override
    protected Integer current() {
      return mX[mIndex];
    }
  }

  public void test() {
    final Iterator<Integer> it = new MockIterator(5, 1, 6, 7, 1, 2, 3, 10);
    assertTrue(it.hasNext());
    assertEquals(Integer.valueOf(6), it.next());
    assertTrue(it.hasNext());
    assertEquals(Integer.valueOf(7), it.next());
    assertTrue(it.hasNext());
    assertEquals(Integer.valueOf(10), it.next());
    assertFalse(it.hasNext());
  }

  public void test1() {
    final Iterator<Integer> it = new MockIterator(5, 6, 7, 1, 2, 3);
    assertTrue(it.hasNext());
    assertEquals(Integer.valueOf(6), it.next());
    assertTrue(it.hasNext());
    assertEquals(Integer.valueOf(7), it.next());
    assertFalse(it.hasNext());
  }

  public void testEmpty() {
    final Iterator<Integer> it = new MockIterator(15, 1, 6, 7, 1, 2, 3, 10);
    assertFalse(it.hasNext());
  }

  public void testEmpty1() {
    final Iterator<Integer> it = new MockIterator(15);
    assertFalse(it.hasNext());
  }
}
