/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
    final ResultStreamHandler rshz = new ResultStreamHandler(mDir, true);
    final File fz = rshz.file("fooblarg");
    final File expz = new File(mDir, "fooblarg");
    assertEquals(expz.getPath(), fz.getPath());
    final ResultStreamHandler rsh = new ResultStreamHandler(mDir, false);
    final File f = rsh.file("fooblarg");
    final File exp = new File(mDir, "fooblarg");
    assertEquals(exp.getPath(), f.getPath());

  }

  /**
   * Test of createFileStream method, of class ResultStreamHandler.
   */
  public void testCreateFileStream() throws IOException {
    final String apples = "apples";
    final ResultStreamHandler rshz = new ResultStreamHandler(mDir, true);
    final File expz = new File(mDir, "foo.gz");
    try (OutputStream fooz = rshz.createFileStream("foo")) {
      fooz.write(apples.getBytes());
    }
    assertEquals(apples, FileHelper.gzFileToString(expz));
    final ResultStreamHandler rsh = new ResultStreamHandler(mDir, false);
    final File exp = new File(mDir, "foo");
    try (OutputStream foo = rsh.createFileStream("foo")) {
      foo.write(apples.getBytes());
    }
    assertEquals(apples, FileUtils.fileToString(exp));
  }


  public void testError() throws IOException {
    assertTrue(mDir.createNewFile());
    final ResultStreamHandler rshz = new ResultStreamHandler(mDir, true);
    try {
      rshz.createFileStream("whatever");
      fail();
    } catch (final IOException e) {

    }

  }
}
