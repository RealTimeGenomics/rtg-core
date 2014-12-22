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
 */
public class CompressedByteArrayTest extends TestCase {

  public void test128() {
    runtestSimple(new CompressedByteArray(30, 128, 1, 7, false));
  }
  public void testDNA() {
    runtestSimple(new CompressedByteArray(30, 5, 3, 7, false));
  }
  public void testProtein() {
    runtestSimple(new CompressedByteArray(30, 22, 2, 9, false));
  }
  public void testBits() {
    runtestSimple(new CompressedByteArray(30, 2, 1, 1, false));
  }

  public void runtestSimple(CompressedByteArray array) {
    final byte[] data = new byte[(int) array.length()];
    for (int i = 0; i < data.length; i++) {
      final int value = (i + 3) % array.getRange();
      assert value < 128;
      data[i] = (byte) value;
    }
    data[data.length - 1] = (byte) (array.getRange() - 1);
    array.set(0, data, data.length);
    final byte[] tmp = new byte[data.length];
    array.get(tmp, 0, data.length);
    for (int i = 0; i < data.length; i++) {
      assertEquals(data[i], tmp[i]);
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
      CompressedByteArray array = new CompressedByteArray(100, 22, 2, 9, false);
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
        }

        // move along to the next set.
        offset += size;
      }
    }
  }

  public void testSetOrder() {
    CompressedByteArray cba = new CompressedByteArray(100, 22, false);
    final byte[] data = {1, 3, 5, 7, 9, 11, 13};
    cba.set(0, data, data.length);
    try {
      cba.set(6, data, data.length);
      fail("expected exception");
    } catch (RuntimeException e) {
      assertEquals("CompressedByteArray.set called out of order", e.getMessage());
    }
  }

  public void testMinBits() {
    assertEquals(1, CompressedByteArray.minBits(0));
    assertEquals(1, CompressedByteArray.minBits(1));
    assertEquals(1, CompressedByteArray.minBits(2));
    assertEquals(2, CompressedByteArray.minBits(4));
    assertEquals(3, CompressedByteArray.minBits(5));
    assertEquals(7, CompressedByteArray.minBits(128));
    assertEquals(31, CompressedByteArray.minBits(Integer.MAX_VALUE));
  }

  public void testDefaultCompression() {
    final int size = 1000;
    // DNA
    CompressedByteArray cba = new CompressedByteArray(size, 5, false);
    long longs = size / 27 + 1;
    assertEquals(longs * 8, cba.bytes());
    assertEquals(size, cba.length());

    // Protein
    cba = new CompressedByteArray(1000, 22, false);
    longs = size / 14 + 1;
    assertEquals(longs * 8, cba.bytes());

    // Quality
    cba = new CompressedByteArray(1000, 64, false);
    longs = size / 10 + 1;
    assertEquals(longs * 8, cba.bytes());

    // 5-bit fields
    cba = new CompressedByteArray(1000, 32, false);
    longs = size / 12 + 1;
    assertEquals(longs * 8, cba.bytes());
  }
}
