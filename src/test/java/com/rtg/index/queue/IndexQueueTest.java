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

import com.rtg.util.array.ExtensibleIndex;
import com.rtg.util.array.longindex.LongChunks;

/**
 */
public class IndexQueueTest extends IndexIteratorQueueTest {

  private static long size(final int bits) {
    return (1L << bits) - 3;
  }

  public void testMemory() {
    checkMemory(200, 0, 8);
    checkMemory(10, 0, 8);
    checkMemory(size(10), 2, 8);
    checkMemory(size(14), 6, 8);
    checkMemory(size(15), 7, 8);
    checkMemory(size(16), 7, 9);
    checkMemory(size(35), 7, 28);
    checkMemory(size(36), 7, 29);
    checkMemory(size(37), 7, 30);
    checkMemory(size(38), 7, 31);
    checkMemory(size(39), 8, 31);
    checkMemory(size(40), 9, 31);
    checkMemory(size(41), 10, 31);
    checkMemory(size(55), 24, 31);
    checkMemory(size(62), 31, 31);
    try {
      checkMemory(-1L, 31, 31);
      fail();
    } catch (final IllegalArgumentException e) {
      // expected
    }
    try {
      checkMemory(Long.MAX_VALUE, 31, 31);
      fail();
    } catch (final IllegalArgumentException e) {
      // expected
    }
  }

  static final int[] X = new int[2];

  void checkMemory(final long length, final int bb, final int cb) {
    new IndexQueue(0, 0, length, 1) {
      @Override
      ExtensibleIndex makeMemory(final int blockBits, final long length, final int subarrayBits) {
        X[0] = blockBits;
        X[1] = subarrayBits;
        return new LongChunks(300); //long enough for all the tests to survive a global integrity check
      }
    };
    //System.err.println(X[0] + " " + X[1]);
    assertEquals(bb, X[0]);
    assertEquals(cb, X[1]);
  }

  public void testMemory2() {
    final ExtensibleIndex makeMemory = new IndexQueue(1, 1, 100, 1).makeMemory(0, 1, 1);
    assertEquals(1, makeMemory.length());
  }

  //case when blocksize might end up being too short
  public void testSmallBlock() {
    final IndexQueue iq = new IndexQueue(1, 7, 100, 1);
    iq.globalIntegrity();
  }

  public void testBug() {
    final IndexQueue iq = new IndexQueue(26, 10, 100, 1);
    final long hash = 129L << 26;
    iq.add(hash, 1337);
    iq.close();
    final QueueIterator qi = iq.iterator(129, 129L << 26);
    assertTrue(qi.next());
    assertEquals(hash, qi.hash());
    assertEquals(1337, qi.id());
    assertFalse(qi.next());
  }
}
