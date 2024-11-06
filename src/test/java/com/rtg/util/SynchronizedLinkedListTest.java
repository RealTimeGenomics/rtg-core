/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
