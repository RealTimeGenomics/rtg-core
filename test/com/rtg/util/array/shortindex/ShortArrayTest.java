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
package com.rtg.util.array.shortindex;


import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;

/**
 * Test Array
 */
public class ShortArrayTest extends AbstractShortIndexTest {

  @Override
  protected ShortIndex create(final long length) {
    return new ShortArray(length);
  }

  @Override
  protected ShortIndex create(final long length, final int bits) {
    //ignore bits
    return new ShortArray(length);
  }


  public void testSerial() throws IOException {
    final ShortArray la = new ShortArray(10);
    for (int i = 0; i < 10; i++) {
      la.set(i, i * 4 + 7);
    }
    final ByteArrayOutputStream out =  new ByteArrayOutputStream();
    la.save(new ObjectOutputStream(out));
    final ByteArrayInputStream in = new ByteArrayInputStream(out.toByteArray());
    final ShortIndex index2 = ShortCreate.loadIndex(new ObjectInputStream(in));
    assertTrue(index2 instanceof ShortArray);
    assertEquals(la.length(), index2.length());
    for (int i = 0; i < 10; i++) {
      assertEquals(la.get(i), index2.get(i));
    }
  }
}

