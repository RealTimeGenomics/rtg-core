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

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.PrintStream;

import com.rtg.util.diagnostic.Diagnostic;

import junit.framework.TestCase;

/**
 * Tests for the corresponding class
 *
 */

public class MultiByteArrayTest extends TestCase {

  public static void main(final String[] args) {
    junit.textui.TestRunner.run(MultiByteArrayTest.class);
  }

  public void test() throws Exception {
    Diagnostic.setLogStream();
    final MultiByteArray mba = new MultiByteArray(50);
    assertEquals(50, mba.length());
    final byte[] b = new byte[3];

    // set 0 bytes
    mba.set(2, new byte[] {3}, 0);
    mba.get(b, 2, 1);
    assertEquals(0, b[0]);

    // set 1 byte
    mba.set(2, new byte[] {3}, 1);
    mba.get(b, 2, 1);
    assertEquals(3, b[0]);

    mba.get(b, 1, 2);
    assertEquals(0, b[0]);
    assertEquals(3, b[1]);
    assertEquals(0, b[2]);
    try {
      mba.load(new ByteArrayInputStream(b), 0, 55);
      fail("expected index out of bounds exception");
    } catch (final IndexOutOfBoundsException ioe) {
      //expected
    }
    mba.load(new ByteArrayInputStream(b), 0, 2);
  }

  public static byte get1(ByteArray bytearray, long offset) {
    final byte[] tmp = new byte[1];
    bytearray.get(tmp, offset, 1);
    return tmp[0];
  }

  public void testMulti() throws IOException {
    final ByteArrayOutputStream bis = new ByteArrayOutputStream();
    Diagnostic.setLogStream(new PrintStream(bis));
    final MultiByteArray mba = new MultiByteArray(18, 3); // 8 bytes per array
    final long len = mba.length();
    assertEquals((long) 18, len);
    assertEquals((byte) 0, get1(mba, len - 1));
    mba.set(len - 1, new byte[] {3}, 1);
    assertEquals(3, get1(mba, len - 1));
    final byte[] data = {1, 4, 9, 16, 25, 36, 49, 64, 81, 100, 121};
    for (int offset = 0; offset < mba.length(); offset += data.length) {
      mba.set(offset, data, Math.min(data.length, (int) (mba.length() - offset)));
    }
    // check byte by byte
    for (long i = mba.length() - 1; i >= 0; i--) {
      assertEquals(data[(int) i % data.length], get1(mba, i));
    }
    // check the multi-getter.
    for (int start = 0; start < 8; start++) {
      byte[] out = new byte[data.length];
      mba.get(out, start, data.length);
      for (int i = 0; i < data.length; i++) {
        assertEquals(data[(i + start) % data.length], out[i]);
      }
    }
    // do a load longer than one segment.
    mba.load(new ByteArrayInputStream(data, 0, data.length), 3L, data.length);
    for (int i = 0; i < data.length; i++) {
      assertEquals(data[i], get1(mba, 3L + i));
    }
    assertTrue(bis.toString().contains("MultiByteArray allocating 8 bytes (block 1 of 3)"));
    assertTrue(bis.toString().contains("MultiByteArray allocating 8 bytes (block 2 of 3)"));
    assertTrue(bis.toString().contains("MultiByteArray allocating 2 bytes (block 3 of 3)"));
  }

  public void testExtend() {
    final MultiByteArray mb = new MultiByteArray(9, 2);
    mb.set(4, (byte) 1);
    mb.set(8, (byte) 2);
    try {
      mb.set(18, (byte) 3);
      fail();
    } catch (ArrayIndexOutOfBoundsException e) {
      //expected
    }
    mb.extendTo(19);
    mb.set(18, (byte) 3);
    final byte[] exp = {0, 0, 0, 0, 1,
                                   0, 0, 0, 2, 0,
                                   0, 0, 0, 0, 0,
                                   0, 0, 0, 3};
    assertEquals((long) exp.length, mb.length());
    for (int i = 0; i < exp.length; i++) {
      assertEquals(exp[i], mb.get(i));
    }
  }
}
