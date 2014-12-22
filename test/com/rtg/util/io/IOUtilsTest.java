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


import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.net.URL;

import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * Tests the IOUtils class. Run from the command line with:<p>
 *
 * java junit.swingui.TestRunner com.reeltwo.util.net.IOUtilsTest java
 * com.reeltwo.util.net.IOUtilsTest
 *
 */
public class IOUtilsTest extends TestCase {

  private static final String STRING = "yollywock";

  private static final byte[] BYTES = STRING.getBytes();

  private static final byte[] EMPTY = new byte[0];


  public void testEmpty() throws IOException {
    testReadAll(EMPTY);
  }


  public void testString() throws IOException {
    testReadAll(BYTES);
  }


  private void checkBytes(final byte[] s, final byte[] bres) {
    assertEquals(s.length, bres.length);
    for (int i = 0; i < s.length; i++) {
      assertEquals(s[i], bres[i]);
    }
  }

  public void testReadAll(final byte[] s) throws IOException {
    InputStream in = new ByteArrayInputStream(s);
    checkBytes(s, IOUtils.readAll(in).getBytes());

    in = new ByteArrayInputStream(s);
    checkBytes(s, IOUtils.readAll(in, "UTF-8").getBytes());
  }

  public void testReadAllFile() throws IOException {
    final File a = FileHelper.createTempFile();
    try {
      PrintStream fw = new PrintStream(a);
      fw.print(STRING);
      fw.close();

      assertEquals(STRING, IOUtils.readAll(a));
      assertEquals(STRING, IOUtils.readAll(new URL(a.toURI().toString())));
      checkBytes(BYTES, IOUtils.readData(new URL(a.toURI().toString())));
    } finally {
      assertTrue(!a.exists() || FileHelper.deleteAll(a));
    }
  }

  public void testReadDataBogusIn() {
    try {
      IOUtils.readData((URL) null);
      fail("Accepted null");
    } catch (final IOException e) {
      fail("IO");
    } catch (final NullPointerException e) {
      // okay
    }

    try {
      IOUtils.readData((InputStream) null);
      fail("Accepted null");
    } catch (final IOException e) {
      fail("IO");
    } catch (final NullPointerException e) {
      // okay
    }
  }


  public void testReadData() {
    String s = "Hobbits live in small holes in the ground";
    try {
      byte[] r = IOUtils.readData(new ByteArrayInputStream(s.getBytes()));
      assertTrue(r != null);
      assertEquals(s, new String(r));
    } catch (final IOException e) {
      fail("IO: " + e.getMessage());
    }
  }


  public void testReadDataEmpty() {
    String s = "";
    try {
      byte[] r = IOUtils.readData(new ByteArrayInputStream(s.getBytes()));
      assertTrue(r != null);
      assertEquals(s, new String(r));
    } catch (final IOException e) {
      fail("IO: " + e.getMessage());
    }
  }


  public void testReadDataZeroLength() {
    try {
      byte[] r = IOUtils.readData(new ByteArrayInputStream(new byte[0]));
      assertTrue(r != null);
      assertEquals(0, r.length);
    } catch (final IOException e) {
      fail("IO: " + e.getMessage());
    }
  }


  public static void main(final String[] args) {
    junit.textui.TestRunner.run(new TestSuite(IOUtilsTest.class));
  }
}
