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

package com.rtg.util;

import java.util.HashSet;
import java.util.LinkedList;

import junit.framework.TestCase;

/**
 */
public class SynchronizedLinkedListTest extends TestCase {

  public void testEverything() {

    final SynchronizedLinkedList<String> sll = new SynchronizedLinkedList<>(new LinkedList<String>());

    assertTrue(sll.add("blah"));
    sll.add("blah2");
    assertFalse(sll.isEmpty());

    final HashSet<String> set = new HashSet<>();
    set.add("blah");
    set.add("blah2");
    assertTrue(sll.containsAll(set));
    assertTrue(sll.contains("blah"));

    assertEquals("blah", sll.removeFirst());
    assertEquals("blah2", sll.remove());

    assertNotNull(sll.toArray());
    assertNotNull(sll.toArray(new String[1]));
    assertNotNull(sll.iterator());
    assertNotNull(sll.descendingIterator());

    assertTrue(sll.addAll(set));

    final HashSet<String> set2 = new HashSet<>();
    set2.add("blah3");
    set2.add("blah4");
    sll.addAll(set2);
    assertEquals(4, sll.size());

    sll.addFirst("blahf");
    assertEquals("blahf", sll.pop());
    sll.addFirst("blahf");
    assertTrue(sll.remove("blahf"));

    assertTrue(sll.retainAll(set));

    assertEquals(2, sll.size());

    assertTrue(sll.offerFirst("blahf"));
    assertTrue(sll.offerLast("blahl"));

    assertEquals("blahf", sll.getFirst());
    assertEquals("blahf", sll.element());
    assertEquals("blahf", sll.removeFirst());
    assertEquals("blahl", sll.getLast());
    assertEquals("blahl", sll.removeLast());
    assertEquals(2, sll.size());

    assertTrue(sll.removeAll(set));
    assertEquals(0, sll.size());

    assertNull(sll.pollFirst());
    assertNull(sll.peekFirst());
    assertNull(sll.pollLast());
    assertNull(sll.peekLast());
    assertTrue(sll.offer("blah"));
    assertTrue(sll.offer("blah"));
    assertTrue(sll.removeFirstOccurrence("blah"));
    assertTrue(sll.removeLastOccurrence("blah"));
    assertFalse(sll.removeLastOccurrence("blah"));
    assertNull(sll.peek());
    assertNull(sll.poll());
  }
}
