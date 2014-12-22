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
package com.rtg.util.bytecompression;

import junit.framework.TestCase;

/**
 * Tests for the corresponding class
 *
 */

public class ByteArrayTest extends TestCase {

  public static void main(final String[] args) {
    junit.textui.TestRunner.run(ByteArrayTest.class);
  }

  protected ByteArray getByteArray(long size, int bits) {
    return ByteArray.allocate(size);
  }

  public void test() {
    final ByteArray ba = ByteArray.allocate(20);
    assertEquals(20, ba.length());
    assertEquals(ba.bytes(), ba.length());
    assertTrue(ba instanceof SingleByteArray);
  }

  public void testSet() {
    final int bits = 3;
    final ByteArray array = getByteArray(128, bits);
    final byte[] data = {1, 3, 5, 7, 0};
    long offset = 0L;
    while (offset + data.length < array.length()) {
      //System.out.println("offset = " + offset);
      array.set(offset, data, data.length);
      assertEquals(1, array.get(offset));
      assertEquals(3, array.get(offset + 1));
      assertEquals(5, array.get(offset + 2));
      assertEquals(7, array.get(offset + 3));
      offset += data.length;
    }
  }

  public void testSimple() {
    final int bits = 3;
    final ByteArray array = getByteArray(128, bits);

    final byte[] data = new byte[(int) array.length()];
    for (int i = 0; i < data.length; i++) {
      final int value = (i + 3) % (1 << bits);
      assert value < 128;
      data[i] = (byte) value;
    }
    data[data.length - 1] = (byte) ((1 << bits) - 1);
    array.set(0L, data, data.length);
    final byte[] tmp = new byte[data.length];
    array.get(tmp, 0L, data.length);
    for (int i = 0; i < data.length; i++) {
      assertEquals("data[" + i + "]", data[i], tmp[i]);
      assertEquals("data[" + i + "]", data[i], array.get(i));
    }
  }

  /**
   * This tests a sequence of sets that moves along the array,
   * and a get that goes back before the set, to make sure nothing
   * was overwritten.
   */
  public void testGetSet() {
    final byte[] data = {1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 20, 18, 16, 14, 12, 10, 8, 6, 4, 2, 0};

    for (int size = 1; size < data.length; size++) {
      ByteArray array = getByteArray(100, 7);
      long offset = 1L;  // because our get goes one byte before offset.
      while (offset < array.length()) {
        final int safeWrite = Math.min(size, (int) (array.length() - offset));
        array.set(offset, data, safeWrite);

        // now check the contents.
        final byte[] tmp = new byte[safeWrite + 1];
        array.get(tmp, offset - 1, safeWrite + 1);
        assertEquals("tmp[0]", offset == 1L ? (byte) 0 : data[size - 1], tmp[0]);
        for (int i = 1; i <= safeWrite; i++) {
          assertEquals("tmp[" + i + "]", data[i - 1], tmp[i]);
          assertEquals("tmp[" + i + "]", data[i - 1], array.get(offset + i - 1));
        }

        // move along to the next set.
        offset += size;
      }
    }
  }

}
