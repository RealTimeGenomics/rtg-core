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
package com.rtg.sam;

import java.util.ArrayList;

import junit.framework.TestCase;


/**
 */
public class CappedConcurrentLinkedListTest extends TestCase {

  public void test() {
    CappedConcurrentLinkedList<Integer> ll = new CappedConcurrentLinkedList<>(1, 0);
    Integer one = 1;
    ll.add(one);
    assertEquals(one, ll.peek());
    assertEquals(one, ll.poll());
    ll.close();
    assertNull(ll.poll());
    assertNull(ll.peek());
    assertEquals("CappedConcurrentLinkedList id: 0 size: 0 hasNext: " + true + " closed: " + true, ll.toString());
    //TODO test all the threading stuff....
  }

  /**
   */
  public void testUnsupported() {
    CappedConcurrentLinkedList<Integer> ll = new CappedConcurrentLinkedList<>(1, 0);
    try {
      ll.addAll(null);
      fail();
    } catch (UnsupportedOperationException e) {
      //expected
    }
    try {
      ll.element();
      fail();
    } catch (UnsupportedOperationException e) {
      //expected
    }
    try {
      ll.iterator();
      fail();
    } catch (UnsupportedOperationException e) {
      //expected
    }
    try {
      ll.offer(null);
      fail();
    } catch (UnsupportedOperationException e) {
      //expected
    }
    try {
      ll.remove();
      fail();
    } catch (UnsupportedOperationException e) {
      //expected
    }
    try {
      ll.remove(null);
      fail();
    } catch (UnsupportedOperationException e) {
      //expected
    }
    try {
      ll.removeAll(null);
      fail();
    } catch (UnsupportedOperationException e) {
      //expected
    }
    try {
      ll.retainAll(null);
      fail();
    } catch (UnsupportedOperationException e) {
      //expected
    }
    //just ignored
    assertFalse(ll.contains(null));
    assertTrue(ll.containsAll(new ArrayList<Integer>()));
    ll.add(1);
    assertEquals(1, ll.size());
    try {
      ll.clear();
      fail();
    } catch (UnsupportedOperationException e) {
      //expected
    }
    assertFalse(ll.isEmpty());
    assertNotNull(ll.toArray());
    assertNotNull(ll.toArray(new Integer[ll.size()]));
  }
}
