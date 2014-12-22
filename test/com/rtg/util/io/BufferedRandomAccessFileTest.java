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

import com.rtg.util.PortableRandom;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * JUnit tests for the BufferedRandomAccessFile class.
 *
 */
public class BufferedRandomAccessFileTest extends TestCase {
  protected File mTmp = null;

  /**
   * Constructor (needed for JUnit)
   *
   * @param name A string which names the object.
   */
  public BufferedRandomAccessFileTest(final String name) {
    super(name);
  }

  @Override
  public void setUp() throws Exception {
    mTmp = File.createTempFile("braf", "test");
  }

  @Override
  public void tearDown() {
    if (mTmp != null) {
      assertTrue(mTmp.delete());
    }
    mTmp = null;

  }

  public void testConstructors() throws Exception {
    try {
      final BufferedRandomAccessFile braf = new BufferedRandomAccessFile(null, "");
      braf.close();
      fail("accepted empty mode");
    } catch (final IllegalArgumentException iae) {
      // expected
    }

    try {
      final BufferedRandomAccessFile braf = new BufferedRandomAccessFile(null, null);
      braf.close();
      fail("accepted null mode");
    } catch (final NullPointerException npe) {
      // expected
    }

    try {
      final BufferedRandomAccessFile braf = new BufferedRandomAccessFile(null, "r");
      braf.close();
      fail("accepted null file");
    } catch (final NullPointerException npe) {
      // expected
    }

    try {
      final BufferedRandomAccessFile braf = new BufferedRandomAccessFile(mTmp, "r");
      braf.close();
    } catch (final Exception e) {
      fail(e.getMessage());
    }

    try {
      final BufferedRandomAccessFile braf = new BufferedRandomAccessFile(mTmp, "r", 0);
      braf.close();
      fail("accepted 0 buffer size");
    } catch (final IllegalArgumentException e) {
      assertTrue(e.getMessage().contains("Buffer size must be >= 1"));
    }

    try {
      final BufferedRandomAccessFile braf = new BufferedRandomAccessFile(mTmp, "r", -10);
      braf.close();
      fail("accepted -ve buffer size");
    } catch (final IllegalArgumentException e) {
      // expected
    }

    try {
      final BufferedRandomAccessFile braf = new BufferedRandomAccessFile(mTmp, "rw");
      braf.close();
      fail("accepted invalid mode");
    } catch (final IllegalArgumentException e) {
      assertTrue(e.getMessage().contains("Can only handle read only mode."));
    }
  }

  private static byte[] setData(final File file) throws IOException {
    final byte[] contents = new byte[1024 * 1027];
    for (int i = 0; i < contents.length; i++) {
      contents[i] = (byte) (i % 256);
    }
    try (RandomAccessFile raf = new RandomAccessFile(file, "rw")) {
      raf.write(contents);
    }
    return contents;
  }

  public void testRead() throws Exception {
    final byte[] contents = setData(mTmp);

    try (RandomAccessFile braf = new BufferedRandomAccessFile(mTmp, "r")) {
      for (int i = 0; i < contents.length; i++) {
        final int x = braf.read();
        assertEquals(contents[i], (byte) x);
        assertEquals(i + 1, braf.getFilePointer());
      }
      assertEquals(-1, braf.read());
      assertEquals(-1, braf.read());

    }
  }

  public void testRead2() throws Exception {
    final byte[] contents = setData(mTmp);

    try (RandomAccessFile braf = new BufferedRandomAccessFile(mTmp, "r")) {
      int count = 0;
      final byte[] buf = new byte[123];
      while (count < contents.length) {
        final int read = braf.read(buf, 0, buf.length);
        assertEquals(count + read, braf.getFilePointer());
        if (read == 0) {
          break;
        }
        for (int i = 0; i < read; i++) {
          assertTrue(i + " : " + count + " : " + read, count + i < contents.length);
          assertEquals("buf[" + i + "]=" + buf[i] + " != contents[" + (count + i) + "]=" + contents[count + i], contents[count + i], buf[i]);
        }
        count += read;
      }
      assertEquals(contents.length, count);
    }
  }

  public void testSeek() throws Exception {
    final byte[] contents = setData(mTmp);

    final byte[] buf = new byte[123];
    int read;

    try (RandomAccessFile braf = new BufferedRandomAccessFile(mTmp, "r")) {
      final PortableRandom rand = new PortableRandom(314159);
      for (int i = 0; i < 100; i++) {
        final int index = rand.nextInt(contents.length - 1000);
        braf.seek(index);
        assertEquals(index, braf.getFilePointer());

        read = braf.read(buf, 0, buf.length);
        assertEquals(buf.length, read);

        for (int j = 0; j < read; j++) {
          assertTrue(j + " : " + index + " : " + read, index + j < contents.length);
          assertEquals("buf[" + j + "]=" + buf[j] + " != contents[" + (index + j) + "]=" + contents[index + j], contents[index + j], buf[j]);
        }
      }

      try {
        braf.seek(-1);
        fail("Accepted negative file offset.");
      } catch (final IOException ioe) {
        // expected
      }

      braf.seek(0);
      assertEquals(0, braf.getFilePointer());
      read = braf.read(buf, 0, buf.length);
      assertEquals(buf.length, read);

      braf.seek(contents.length + 10);
      assertEquals(contents.length + 10, braf.getFilePointer());
      read = braf.read(buf, 0, buf.length);
      assertEquals(-1, read);

      braf.seek(2 * contents.length + 10);
      assertEquals(2 * contents.length + 10, braf.getFilePointer());
      read = braf.read(buf, 0, buf.length);
      assertEquals(-1, read);

      read = braf.read(buf, 0, buf.length);
      assertEquals(-1, read);

    }
  }

  public void testSeekSame() throws Exception {
    final byte[] contents = setData(mTmp);

    final byte[] buf = new byte[123];
    int read;

    try (RandomAccessFile braf = new BufferedRandomAccessFile(mTmp, "r")) {
      braf.seek(0);
      //assertEquals(0, braf.getFilePointer());
      read = braf.read(buf, 0, buf.length);
      assertEquals(buf.length, read);
      for (int i = 0; i < buf.length; i++) {
        assertEquals(contents[i], buf[i]);
      }

      braf.seek(0);
      //assertEquals(0, braf.getFilePointer());
      read = braf.read(buf, 0, buf.length);
      assertEquals(buf.length, read);
      for (int i = 0; i < buf.length; i++) {
        assertEquals(contents[i], buf[i]);
      }

      braf.seek(10);
      assertEquals(10, braf.getFilePointer());
      read = braf.read(buf, 0, buf.length);
      assertEquals(buf.length, read);
      for (int i = 0; i < buf.length; i++) {
        assertEquals(contents[i + 10], buf[i]);
      }

      braf.seek(0);
      //assertEquals(0, braf.getFilePointer());
      read = braf.read(buf, 0, buf.length);
      assertEquals(buf.length, read);
      for (int i = 0; i < buf.length; i++) {
        assertEquals(contents[i], buf[i]);
      }

      braf.seek(100000);
      assertEquals(100000, braf.getFilePointer());
      read = braf.read(buf, 0, buf.length);
      assertEquals(buf.length, read);
      for (int i = 0; i < buf.length; i++) {
        assertEquals(contents[i + 100000], buf[i]);
      }

      braf.seek(100001);
      assertEquals(100001, braf.getFilePointer());
      read = braf.read(buf, 0, buf.length);
      assertEquals(buf.length, read);
      for (int i = 0; i < buf.length; i++) {
        assertEquals(contents[i + 100001], buf[i]);
      }

      braf.seek(99999);
      assertEquals(99999, braf.getFilePointer());
      read = braf.read(buf, 0, buf.length);
      assertEquals(buf.length, read);
      for (int i = 0; i < buf.length; i++) {
        assertEquals(contents[i + 99999], buf[i]);
      }
    }
  }

  public static Test suite() {
    return new TestSuite(BufferedRandomAccessFileTest.class);
  }

  /**
   * Main to run from tests from command line.
   *
   * @param args ignored.
   */
  public static void main(final String[] args) {
    junit.textui.TestRunner.run(suite());
  }

}
