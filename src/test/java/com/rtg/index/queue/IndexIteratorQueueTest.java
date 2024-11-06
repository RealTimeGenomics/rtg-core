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
