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

import com.rtg.util.array.ExtensibleIndex;
import com.rtg.util.array.longindex.LongChunks;



/**
 */
public class IndexQueueTest extends IndexIteratorQueueTest {

  private static long size(final int bits) {
    return (1L << bits) - 2;
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
      checkMemory(Long.MAX_VALUE, 31, 31);
      fail();
    } catch (final RuntimeException e) {
      // expected
      assertEquals("Length too large:" + Long.MAX_VALUE + " for radix bits=0", e.getMessage());
    }
    try {
      checkMemory(Long.MAX_VALUE - 2, 31, 31);
      fail();
    } catch (final IllegalArgumentException e) {
      // expected
      assertEquals("Number out of range:" + (Long.MAX_VALUE - 1), e.getMessage());
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
