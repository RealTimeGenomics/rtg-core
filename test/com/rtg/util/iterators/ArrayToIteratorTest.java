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
public class ArrayToIteratorTest extends TestCase {

  private static final MockTransform MOCK_TRANSFORM = new MockTransform();

  private static class MockTransform extends Transform<String, Integer> {
    @Override
    public Integer trans(String x) {
      return Integer.valueOf(x);
    }
  }

  public void testEmpty() {
    final Iterator<String> it = Transform.array2Iterator(new String[0]);
    final Iterator<Integer> trans = MOCK_TRANSFORM.trans(it);
    final Iterator<Integer> exp = Transform.array2Iterator(new Integer[0]);
    TransformTest.check(exp , trans);
  }

  public void testArray2Iterator() {
    final Iterator<String> act = Transform.array2Iterator(new String[] {"1", "3", "2"});
    assertTrue(act.hasNext());
    assertEquals("1", act.next());
    assertTrue(act.hasNext());
    assertEquals("3", act.next());
    assertTrue(act.hasNext());
    assertEquals("2", act.next());
    assertFalse(act.hasNext());
  }
}
