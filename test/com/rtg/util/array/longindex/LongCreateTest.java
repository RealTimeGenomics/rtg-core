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
package com.rtg.util.array.longindex;


import junit.framework.TestCase;

/**
 * Test Create
 */
public class LongCreateTest extends TestCase {

  private static final long FREE_LIMIT = 8L * Integer.MAX_VALUE + 1000000000L; //allow 1000m of freeboard

  public void testBad() {
    try {
      LongCreate.createIndex(-1L);
      fail("NegativeArraySizeException expected");
    } catch (final NegativeArraySizeException e) {
      //expected
      assertEquals("Negative length=-1", e.getMessage());
    }

    final LongIndex index = LongCreate.createIndex(0);
    index.integrity();
  }

  public void test() {
    final LongIndex a = LongCreate.createIndex(10);
    assertEquals(10, a.length());
    assertTrue(a instanceof LongArray);

    //Only run if there is enough memory not to nuke everything
    System.gc();
    final long mem = Runtime.getRuntime().freeMemory();
    //System.err.println("FREE_LIMIT=" + FREE_LIMIT + " mem=" + mem);
    if (mem > FREE_LIMIT && Runtime.getRuntime().freeMemory() < 4000000000L) {
      // SAI: There can be some JVM overheads that prevent Integer.MAX_VALUE lengths
      final int safeLength = Integer.MAX_VALUE - 5;
      final LongIndex b = LongCreate.createIndex(safeLength);
      assertEquals(safeLength, b.length());
      assertTrue(b instanceof LongArray);
      System.gc();

      final LongIndex c = LongCreate.createIndex(Integer.MAX_VALUE + 1L);
      assertEquals(Integer.MAX_VALUE + 1L, c.length());
      assertTrue(c instanceof LongChunks);
      System.gc();
    }
  }
}
