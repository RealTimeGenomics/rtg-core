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

import java.util.Iterator;

import junit.framework.TestCase;

/**
 */
public class BasicLinkedListNodeTest extends TestCase {

  private BasicLinkedListNode<String> getList(String[] values) {
    BasicLinkedListNode<String> head = null;
    for (int i = values.length - 1; i >= 0; i--) {
      head = new BasicLinkedListNode<>(values[i], head);
    }
    return head;
  }

  public void test() {
    String[] values = {"a", "b", "c", "d"};
    BasicLinkedListNode<String> f = getList(values);
    assertEquals("a", f.getValue());
    assertEquals(4, f.size());
    Iterator<String> it = f.iterator();
    for (String s : values) {
      assertTrue(it.hasNext());
      assertEquals(s, it.next());
    }
    assertFalse(it.hasNext());

    BasicLinkedListNode<String> firstSide = new BasicLinkedListNode<>("first", f);
    BasicLinkedListNode<String> secondSide = new BasicLinkedListNode<>("second", f);
    assertEquals(5, firstSide.size());
    assertEquals(5, secondSide.size());

    it = firstSide.iterator();
    assertTrue(it.hasNext());
    assertEquals("first", it.next());
    for (String s : values) {
      assertTrue(it.hasNext());
      assertEquals(s, it.next());
    }
    assertFalse(it.hasNext());

    it = secondSide.iterator();
    assertTrue(it.hasNext());
    assertEquals("second", it.next());
    for (String s : values) {
      assertTrue(it.hasNext());
      assertEquals(s, it.next());
    }
    assertFalse(it.hasNext());
  }

  public void testNotSupported() {
    BasicLinkedListNode<String> list = new BasicLinkedListNode<>("second", null);
    Iterator<String> it = list.iterator();
    try {
      it.remove();
      fail("Expected exception not thrown");
    } catch (UnsupportedOperationException e) {
      assertEquals("Remove not supported", e.getMessage());
    }
  }

}
