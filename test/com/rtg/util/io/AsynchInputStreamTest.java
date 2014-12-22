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
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;

import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;


/**
 */
public class AsynchInputStreamTest extends TestCase {

  private static final String EXAMPLE1 = "abc\ndef ghi\n";

  AsynchInputStream getStream(File file, String text) throws IOException {
    if (text != null) {
      FileUtils.stringToFile(text, file);
    }
    return new AsynchInputStream(new FileInputStream(file));
  }

  public void testReadSmall() throws IOException {
    File file = File.createTempFile("test", "gzipasynch");
    try {
      try (AsynchInputStream input = getStream(file, EXAMPLE1)) {
        final byte[] buf = new byte[100];
        assertEquals(12, input.read(buf));
        for (int i = 0; i < EXAMPLE1.length(); i++) {
          assertEquals(EXAMPLE1.charAt(i), (char) buf[i]);
        }
        final int exp = -1;
        assertEquals(exp, input.read(buf));
      }
    } finally {
      assertTrue(FileHelper.deleteAll(file));
    }
  }

  public void testEarlyClose() throws IOException {
    Diagnostic.setLogStream();
    File file = File.createTempFile("test", "gzipasynch");
    try {
      AsynchInputStream input = getStream(file, EXAMPLE1);
      input.close(); // check that this does not throw an exception
      // Note: it would be nice to check that the close has stopped the
      // input from being read, but we do not know how far the thread has
      // already read, so there is a race condition if we uncomment the
      // following lines.  (It fails on some linux machines, but not others).
      // Is there any nice way of checking that the early close is working?
      // final byte[] buf = new byte[100];
      // assertEquals(-1, input.read(buf));
    } finally {
      assertTrue(FileHelper.deleteAll(file));
    }
  }

  public void testMarkNotSupported() throws IOException {
    Diagnostic.setLogStream();
    File file = File.createTempFile("test", "gzipasynch");
    try {
      try (AsynchInputStream input = getStream(file, EXAMPLE1)) {
        assertFalse(input.markSupported());
      }
    } finally {
      assertTrue(FileHelper.deleteAll(file));
    }
  }

  public void testEmpty() throws IOException {
    Diagnostic.setLogStream();
    final File file = File.createTempFile("test", "gzipasynch");
    try {
      try (AsynchInputStream in = getStream(file, null)) {
        assertEquals(1024 * 1024, in.mQueue.maxSize());
        final byte[] buf = new byte[1];
        assertEquals(-1, in.read(buf, 0, 1));
      }
    } finally {
      assertTrue(FileHelper.deleteAll(file));
    }
  }

  //Testing the Self-suppression problem
  public void testExceptionHandling() throws InterruptedException {
    final byte[] buff = new byte[3];
    try {
      try (final AsynchInputStream stream = new AsynchInputStream(new InputStream() {
        int mNum = 2;
        @Override
        public int read() throws IOException {
          if (mNum > 0) {
            return mNum--;
          } else {
            throw new IOException("Expected");
          }
        }
      })) {
        final int b = stream.read(buff, 0, 3);
        final int n = stream.read(buff, 0, 1);
        //These assert statements should not be reached, the second read statement should be throwing an exception
        //They are here because findbugs does not like ignoring the return value of a method
        assertEquals(2, b);
        assertEquals(-1, n);
      }
      fail("Should have thrown an IOException");
    } catch (IOException e) {
      assertEquals("Expected", e.getMessage());
    }
  }

}
