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
package com.rtg.util;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;

import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;

/**
 * Test class
 */
public class ResultStreamHandlerTest extends TestCase {

  private File mDir = null;

  @Override
  public void setUp() throws IOException {
    mDir = FileUtils.createTempDir("resultstreamhandler", "test");
    assertTrue(mDir.delete());
  }

  @Override
  public void tearDown() {
    assertTrue(mDir == null || !mDir.exists() || FileHelper.deleteAll(mDir));
    mDir = null;
  }

  /**
   * Test of file method, of class ResultStreamHandler.
   */
  public void testFile() {
    final ResultStreamHandler rshz = new ResultStreamHandler(mDir, true, ".gzip");
    final File fz = rshz.file("fooblarg");
    final File expz = new File(mDir, "fooblarg");
    assertEquals(expz.getPath(), fz.getPath());
    final ResultStreamHandler rsh = new ResultStreamHandler(mDir, false, ".gzip");
    final File f = rsh.file("fooblarg");
    final File exp = new File(mDir, "fooblarg");
    assertEquals(exp.getPath(), f.getPath());

  }

  /**
   * Test of createFileStream method, of class ResultStreamHandler.
   */
  public void testCreateFileStream() throws IOException {
    final String apples = "apples";
    final ResultStreamHandler rshz = new ResultStreamHandler(mDir, true, ".gzip");
    final File expz = new File(mDir, "foo.gzip");
    try (OutputStream fooz = rshz.createFileStream("foo")) {
      fooz.write(apples.getBytes());
    }
    assertEquals(apples, FileHelper.gzFileToString(expz));
    final ResultStreamHandler rsh = new ResultStreamHandler(mDir, false, ".gzip");
    final File exp = new File(mDir, "foo");
    try (OutputStream foo = rsh.createFileStream("foo")) {
      foo.write(apples.getBytes());
    }
    assertEquals(apples, FileUtils.fileToString(exp));
  }


  public void testError() throws IOException {
    assertTrue(mDir.createNewFile());
    final ResultStreamHandler rshz = new ResultStreamHandler(mDir, true, ".gzip");
    try {
      rshz.createFileStream("whatever");
      fail();
    } catch (final IOException e) {

    }

  }
}
