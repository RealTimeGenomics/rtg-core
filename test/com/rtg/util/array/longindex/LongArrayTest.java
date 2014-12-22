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

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;

/**
 * Test Array
 */
public class LongArrayTest extends AbstractLongIndexTest {

  @Override
  protected LongIndex create(final long length) {
    return new LongArray(length);
  }

  @Override
  protected LongIndex create(final long length, final int bits) {
    //ignore bits
    return new LongArray(length);
  }

  public void testBadLengthExtra() {
    final LongArray a = new LongArray(5);
    final long value = (long) Integer.MAX_VALUE + 2L;
    try {
      a.get(value);
      fail();
    } catch (final IndexOutOfBoundsException ioobe) {
      assertEquals("ii=-2147483647 mArrays.length=5 index=2147483649", ioobe.getMessage());
    }
  }

  public void testSerial() throws IOException {
    final LongArray la = new LongArray(10);
    for (int i = 0; i < 10; i++) {
      la.set(i, i * 4 + 7);
    }
    final ByteArrayOutputStream out =  new ByteArrayOutputStream();
    la.save(new ObjectOutputStream(out));
    final ByteArrayInputStream in = new ByteArrayInputStream(out.toByteArray());
    final LongIndex index2 = LongCreate.loadIndex(new ObjectInputStream(in));
    assertTrue(index2 instanceof LongArray);
    assertEquals(la.length(), index2.length());
    for (int i = 0; i < 10; i++) {
      assertEquals(la.get(i), index2.get(i));
    }
  }
}

