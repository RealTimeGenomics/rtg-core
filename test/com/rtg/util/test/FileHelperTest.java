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
package com.rtg.util.test;

import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.IOException;
import java.io.StringReader;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * Tests for FileHelper.
 *
 */
public class FileHelperTest extends TestCase {

  public FileHelperTest(final String name) {
    super(name);
  }

  public void testCreateTempFile() throws Exception {
    final File f = FileHelper.createTempFile();
    try {
      assertNotNull(f);
      assertTrue(f.exists());
      assertTrue(f.isFile());
      assertTrue(f.canRead());
      assertEquals(0, f.length());
      assertTrue(f.getName().startsWith("unit"));
      assertTrue(f.getName().endsWith("test"));
    } finally {
      assertTrue(f.delete());
    }
    assertFalse(f.exists());

    final File f2 = FileHelper.createTempFile();
    try {
      assertNotNull(f2);
      assertTrue(f2.exists());
      assertTrue(f2.isFile());
      assertTrue(f2.canRead());
      assertEquals(0, f2.length());

      assertTrue(!f2.equals(f));
      assertTrue(f2.getName().startsWith("unit"));
      assertTrue(f2.getName().endsWith("test"));
    } finally {
      assertTrue(f2.delete());
    }
    assertFalse(f2.exists());
  }

  public void testCreateTempDirectory() throws Exception {
    final File f = FileHelper.createTempDirectory();
    try {
      assertNotNull(f);
      assertTrue(f.exists());
      assertTrue(f.isDirectory());
      assertTrue(f.canRead());
      assertTrue(f.getName().startsWith("unit"));
      assertTrue(f.getName().endsWith("test"));
    } finally {
      assertTrue(FileHelper.deleteAll(f));
    }
    assertFalse(f.exists());

    final File f2 = FileHelper.createTempDirectory();
    try {
      assertNotNull(f2);
      assertTrue(f2.exists());
      assertTrue(f2.isDirectory());
      assertTrue(f2.canRead());

      assertTrue(!f2.equals(f));
      assertTrue(f2.getName().startsWith("unit"));
      assertTrue(f2.getName().endsWith("test"));
    } finally {
      assertTrue(FileHelper.deleteAll(f2));
    }
    assertFalse(f2.exists());
  }

  public void testSubDirMethods() throws IOException {
    final File parent = FileHelper.createTempDirectory();
    try {
      final File f = FileHelper.createTempFile(parent);
      assertTrue(f.getName().startsWith("unit"));
      assertTrue(f.getName().endsWith("test"));
      assertEquals(f.getParentFile(), parent);
      final File f2 = FileHelper.createTempDirectory(parent);
      assertTrue(f2.getName().startsWith("unit"));
      assertTrue(f2.getName().endsWith("test"));
      assertEquals(f2.getParentFile(), parent);
    } finally {
      assertTrue(FileHelper.deleteAll(parent));
    }
  }

  public void testReaderToString() throws IOException {
    assertEquals("", FileHelper.readerToString(new StringReader("")));
    assertEquals(" ", FileHelper.readerToString(new StringReader(" ")));
  }

  public void testResourceToString() throws IOException {
    final String str = FileHelper.resourceToString("com/rtg/util/FileUtilsCartesianTest.properties");
    assertTrue(str.contains("this as a resource"));
  }

  public void testResourceToStringBad() throws IOException {
    try {
      FileHelper.resourceToString("bad");
      fail();
    } catch (final RuntimeException e) {
      // expected
      assertEquals("Unable to find resource:bad", e.getMessage());
    }
  }

  public void testStreamToFile() throws IOException {
    final File f = FileHelper.createTempFile();
    try {
      assertEquals(f, FileHelper.streamToFile(new ByteArrayInputStream(new byte[0]), f));
      assertEquals(0, f.length());
      assertEquals(f, FileHelper.streamToFile(new ByteArrayInputStream(new byte[42]), f));
      assertEquals(42, f.length());
      try {
        FileHelper.streamToFile(null, f);
        fail();
      } catch (final NullPointerException e) {
        assertEquals("null stream given", e.getMessage());
      }
      try {
        FileHelper.streamToFile(new ByteArrayInputStream(new byte[0]), null);
        fail();
      } catch (final NullPointerException e) {
        assertEquals("null file given", e.getMessage());
      }
    } finally {
      assertTrue(FileHelper.deleteAll(f));
    }
  }

  public void testStreamToGzFile() throws IOException {
    final File f = FileHelper.createTempFile();
    try {
      assertEquals(f, FileHelper.streamToGzFile(new ByteArrayInputStream(new byte[0]), f));
      assertTrue(f.exists());
      assertEquals(f, FileHelper.streamToGzFile(new ByteArrayInputStream(new byte[42]), f));
      assertTrue(f.length() > 0);
      assertEquals(42, FileHelper.gzFileToString(f).length());
      try {
        FileHelper.streamToFile(null, f);
        fail();
      } catch (final NullPointerException e) {
        assertEquals("null stream given", e.getMessage());
      }
      try {
        FileHelper.streamToFile(new ByteArrayInputStream(new byte[0]), null);
        fail();
      } catch (final NullPointerException e) {
        assertEquals("null file given", e.getMessage());
      }
    } finally {
      assertTrue(FileHelper.deleteAll(f));
    }
  }

  public void testStringToGzFile() throws IOException {
    final File f = FileHelper.createTempFile();
    try {
      assertEquals(f, FileHelper.stringToGzFile("", f));
      assertEquals("", FileHelper.gzFileToString(f));
      assertEquals(f, FileHelper.stringToGzFile("hi", f));
      assertEquals("hi", FileHelper.gzFileToString(f));
      try {
        FileHelper.stringToGzFile(null, f);
        fail();
      } catch (final NullPointerException e) {
        assertEquals("null string given", e.getMessage());
      }
    } finally {
      assertTrue(FileHelper.deleteAll(f));
    }
  }

  public void testDeleteAllBad() {
    assertTrue(FileHelper.deleteAll(new File("there-is-no-file-called-this-i-hope")));
  }

  /**
   * Adds tests to suite to be run by main
   * @return The test suite.
   */
  public static Test suite() {
    return new TestSuite(FileHelperTest.class);
  }

  /**
   * Main to run from tests from command line.
   * @param args ignored.
   */
  public static void main(final String[] args) {
    junit.textui.TestRunner.run(suite());
  }


}
