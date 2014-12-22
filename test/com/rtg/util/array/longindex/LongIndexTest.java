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

import java.io.IOException;
import java.io.ObjectOutputStream;

/**
 * Test Chunks
 */
public class LongIndexTest extends LongChunksTest {

  @Override
  protected LongIndex create(final long length) {
    return new LongChunks(length, 5);
  }

  public void testStuff() {

    try {
      new LongIndex(-1L) {

        @Override
        public void set(final long index, final long value) { }

        @Override
        public long get(final long index) {
          return 0;
        }

        @Override
        public boolean safeFromWordTearing() {
          return true;
        }

        @Override
        public void save(ObjectOutputStream dos) throws IOException {
          throw new UnsupportedOperationException("Not implemented yet");
        }
      };
      fail();
    } catch (final NegativeArraySizeException nase) {
      assertEquals("length=-1", nase.getMessage());
    }

    final LongIndex li = new LongIndex(4) {
      @Override
      public void set(final long index, final long value) { }
      @Override
      public long get(final long index) {
        return 5;
      }

      @Override
      public boolean safeFromWordTearing() {
        return true;
      }

      @Override
      public void save(ObjectOutputStream dos) throws IOException {
        throw new UnsupportedOperationException("Not implemented yet");
      }
    };
    li.set(1, 5);
    assertEquals(5L, li.get(1));
    assert li.integrity();
  }
}
