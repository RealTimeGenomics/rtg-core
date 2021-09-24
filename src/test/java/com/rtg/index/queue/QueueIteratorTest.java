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

//import com.rtg.util.array.longindex.LongChunks;

/**
 */
public class QueueIteratorTest extends IndexIteratorQueueTest {
  public void testIt1() {
    //final LongChunks lc = new LongChunks(4, 100, 5);
    final IndexQueue iq = new IndexQueue(3, 2, 100, 2);
    iq.add(2, 3);
    iq.close();
    final QueueIterator it = iq.iterator(0, 0);
    it.integrity();
    assertEquals("QueueIterator:", it.toString());
    assertTrue(it.next());
    assertEquals(2, it.hash());
    assertEquals(3, it.id());
    assertFalse(it.next());
    try {
      it.hash();
      fail();
    } catch (final IllegalStateException e) {
      // expected
    }
    try {
      it.id();
      fail();
    } catch (final IllegalStateException e) {
      // expected
    }

    assertFalse(iq.iterator(1, 1 << 3).next());
  }
}
