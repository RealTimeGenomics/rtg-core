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
package com.rtg.util.array.objectindex;


import junit.framework.Test;
import junit.framework.TestSuite;

/**
 * Test Chunks
 */
public class ObjectChunksTest extends AbstractObjectIndexTest {

  private static final int CHUNK_BITS = 29;

  public static Test suite() {
    return new TestSuite(ObjectChunksTest.class);
  }

  /**
   * Constructor for ChunksTest.
   */
  public ObjectChunksTest(final String arg0) {
    super(arg0);
  }

  @Override
  protected ObjectIndex<Integer> create(final long length) {
    return new ObjectChunks<>(length);
  }

  @Override
  protected ObjectIndex<Integer> create(final long length, final int bits) {
    return new ObjectChunks<>(length, bits);
  }

  public void testChunkSize() {
    final ObjectChunks<Integer> dc = new ObjectChunks<>(100L);
    dc.integrity();
    assertEquals(100, dc.length());
    assertEquals(1L << 28, dc.chunkSize());
  }

  public void testTooLong() {
    if (Runtime.getRuntime().freeMemory() < 4000000000L) {
      try {
        new ObjectChunks<Integer>((1L << CHUNK_BITS) * Integer.MAX_VALUE + 1L, CHUNK_BITS);
        fail("RuntimeException expected");
      } catch (final RuntimeException e) {
        //expected
      }
      try {
        new ObjectChunks<Integer>((1L << CHUNK_BITS - 1) * Integer.MAX_VALUE + 1L, CHUNK_BITS);
        fail("OutOfMemoryError expected");
      } catch (final OutOfMemoryError e) {
        //expected
      }
    }
  }

  /** Tests allocations of chunks and formula for computing how many chunks there are. */
  public void testEdgeCase() {
      final int bits = 3;
      final int size = 1 << bits;
      final ObjectIndex<Integer> ix0 = new ObjectChunks<>(size - 1, bits);
      ix0.integrity();
      final ObjectIndex<Integer> ix1 = new ObjectChunks<>(size, bits);
      ix1.integrity();
      final ObjectIndex<Integer> ix2 = new ObjectChunks<>(size + 1, bits);
      ix2.integrity();
  }
}
