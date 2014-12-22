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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

import junit.framework.TestCase;

/**
 */
public class ComposeIteratorsTest extends TestCase {

  private static class Trans extends Transform<List<Integer>, Iterator<Integer>> {
    @Override
    public Iterator<Integer> trans(List<Integer> x) {
      return x.iterator();
    }
  }

  public void test() {
    final ArrayList<List<Integer>> lists = new ArrayList<>();
    lists.add(Arrays.asList(0, 1, 2, 3, 4));
    lists.add(Arrays.asList(5));
    lists.add(new ArrayList<Integer>());
    lists.add(Arrays.asList(6, 7));
    //final Iterator<Integer> iterator = new ComposeIterators<>(lists.iterator(), trans);
    final Iterator<Integer> iterator = Transform.flatten(lists.iterator());
    for (int i = 0; i < 8; i++) {
      assertTrue("" + i, iterator.hasNext());
      assertEquals(i, iterator.next().intValue());

    }
    assertFalse(iterator.hasNext());
  }

  public void testFlatten2() {
    final ArrayList<List<Integer>> lists = new ArrayList<>();
    lists.add(Arrays.asList(0, 1, 2, 3, 4));
    lists.add(Arrays.asList(5));
    lists.add(new ArrayList<Integer>());
    lists.add(Arrays.asList(6, 7));
    //final Iterator<Integer> iterator = new ComposeIterators<>(lists.iterator(), trans);
    final Iterator<Integer> iterator = Transform.flatten(lists.iterator(), new Trans());
    for (int i = 0; i < 8; i++) {
      assertTrue("" + i, iterator.hasNext());
      assertEquals(i, iterator.next().intValue());

    }
    assertFalse(iterator.hasNext());
  }

  public void testEmpty() {
    final Iterator<Integer> iterator = Transform.flatten(new ArrayList<List<Integer>>().iterator(), new Trans());
    assertFalse(iterator.hasNext());
  }

  public void testEmptyInternals() {
    final ArrayList<List<Integer>> lists = new ArrayList<>();
    lists.add(new ArrayList<Integer>());
    lists.add(new ArrayList<Integer>());
    final Iterator<Integer> iterator = Transform.flatten(lists.iterator(), new Trans());
    assertFalse(iterator.hasNext());
  }

  public void testEmptyFirst() {
    final ArrayList<List<Integer>> lists = new ArrayList<>();
    lists.add(new ArrayList<Integer>());
    lists.add(Arrays.asList(0, 1, 2));
    final Iterator<Integer> iterator = Transform.flatten(lists.iterator(), new Trans());
    for (int i = 0; i < 3; i++) {
      assertTrue("" + i, iterator.hasNext());
      assertEquals(i, iterator.next().intValue());

    }
    assertFalse(iterator.hasNext());
  }
}
