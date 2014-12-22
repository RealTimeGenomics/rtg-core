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

package com.rtg.util.io;

import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.Arrays;

import junit.framework.TestCase;

/**
 * Test class
 */
public class RandomAccessFileStreamTest extends TestCase {

  public void testSomeMethod() throws IOException {
    final File tmp = FileUtils.stringToFile("01234567890123456789", File.createTempFile("tmp", "file"));
    try {
      try (RandomAccessFileStream raf = new RandomAccessFileStream(new RandomAccessFile(tmp, "r"))) {
        raf.seek(5);
        assertEquals('5', (char) raf.read());
        raf.seek(11);
        final byte[] buf = new byte[5];
        final int len = raf.read(buf);
        assertEquals("12345".substring(0, len), new String(buf, 0, len));

        assertEquals(16, raf.getPosition());

        assertEquals(20, raf.length());

        assertEquals(4, raf.available());

        assertFalse(raf.markSupported());

        final byte[] b = new byte[2];
        assertEquals(2, raf.read(b));
        assertTrue(Arrays.equals(new byte[]{(byte) '6', (byte) '7'}, b));
        assertEquals(2, raf.read(b, 0, 2));
        assertTrue(Arrays.toString(b), Arrays.equals(new byte[]{(byte) '8', (byte) '9'}, b));
      }
    } finally {
      assertTrue(tmp.delete());
    }
  }
}
