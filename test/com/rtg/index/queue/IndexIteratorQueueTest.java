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
package com.rtg.index.queue;

import junit.framework.TestCase;


/**
 */
public class IndexIteratorQueueTest extends TestCase {

  public void test1() {
    final IndexQueue iq = new IndexQueue(3, 2, 100, 2);
    assertEquals("IndexQueue lower=3 upper=2", iq.toString());
    iq.globalIntegrity();
    iq.add(2, 3);
    iq.globalIntegrity();
    iq.close();
    try {
      iq.add(2, 3);
      fail();
    } catch (final IllegalStateException e) {
      //expected
    }
    iq.globalIntegrity();
    final QueueIterator it = iq.iterator(0, 0);
    assertTrue(it.next());
    assertEquals(2, it.hash());
    assertEquals(3, it.id());

    assertFalse(iq.iterator(1, 1 << 3).next());
  }

  public void test2() {
    final IndexQueue iq = new IndexQueue(3, 2, 100, 8);
    iq.globalIntegrity();
    for (int i = 0; i < 200; ++i) {
      iq.add(8, i + 1);
    }
    iq.globalIntegrity();
    iq.close();
    iq.globalIntegrity();
    final QueueIterator it = iq.iterator(1, 1 << 3);
    for (int i = 0; i < 200; ++i) {
      assertTrue(it.next());
      assertEquals(8, it.hash());
      assertEquals(i + 1, it.id());
    }
    assertFalse(it.next());
  }

  public void test3() {
    final IndexQueue iq = new IndexQueue(3, 2, 100, 5);
    iq.globalIntegrity();
    for (int r = 0; r < 4; ++r) {
      for (int i = 0; i < 30; ++i) {
        iq.add(r << 3, i + 1);
      }
    }
    iq.globalIntegrity();
    iq.close();
    iq.globalIntegrity();
    for (int r = 0; r < 4; ++r) {
      final QueueIterator it = iq.iterator(r, r << 3);
      for (int i = 0; i < 30; ++i) {
        assertTrue(it.next());
        assertEquals(r << 3, it.hash());
        assertEquals(i + 1, it.id());
      }
      assertFalse(it.next());
    }
  }

}
