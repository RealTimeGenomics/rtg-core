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
import java.util.NoSuchElementException;

import junit.framework.TestCase;

/**
 */
public class TransformTest extends TestCase {

  private static final MockTransform MOCK_TRANSFORM = new MockTransform();

  private static final MockInt MOCK_INT = new MockInt();

  private static class MockTransform extends Transform<String, Integer> {
    @Override
    public Integer trans(String x) {
      return Integer.valueOf(x);
    }
  }

  private static class MockInt extends Transform<Integer, Integer> {
    @Override
    public Integer trans(Integer x) {
      return -x;
    }
  }

  static <X> void check(Iterator<X> exp, Iterator<X> actual) {
    while (exp.hasNext()) {
      assertTrue(actual.hasNext());
      assertEquals(exp.next(), actual.next());
    }
    assertFalse(actual.hasNext());
  }

  public void test() {
    final Iterator<String> it = Transform.array2Iterator(new String[] {"1", "3", "2"});
    final Iterator<Integer> trans = MOCK_TRANSFORM.trans(it);
    final Iterator<Integer> exp = Transform.array2Iterator(new Integer[] {1, 3, 2});
    check(exp , trans);
  }

  public void testArray2IteratorEmpty() {
    final Iterator<String> act = Transform.array2Iterator(new String[0]);
    assertFalse(act.hasNext());
    try {
      act.next();
      fail();
    } catch (final NoSuchElementException e) {
      //expected
    }
  }

  public void testCompose() {
    final Transform<String, Integer> comp = Transform.compose(MOCK_TRANSFORM, MOCK_INT);
    assertEquals(Integer.valueOf(-1), comp.trans("1"));
  }
}
