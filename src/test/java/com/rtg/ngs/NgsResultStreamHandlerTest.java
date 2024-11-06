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
package com.rtg.ngs;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;

import com.rtg.AbstractTest;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;

/**
 * Test class
 */
public class NgsResultStreamHandlerTest extends AbstractTest {

  File mDir = null;

  @Override
  public void setUp() throws IOException {
    super.setUp();
    mDir = FileUtils.createTempDir("ngsresultstreamhandler", "test");
    assertTrue(mDir.delete());
  }

  @Override
  public void tearDown() throws IOException {
    assertTrue(mDir == null || !mDir.exists() || FileHelper.deleteAll(mDir));
    mDir = null;
    super.tearDown();
  }

  private void checkStream(OutputStream out, String name, boolean zip) throws IOException {
    final String apples = "apples";
    try {
      out.write(apples.getBytes());
    } finally {
      out.close();
    }
    final File f = new File(mDir, name + (zip ? FileUtils.GZ_SUFFIX : ""));
    final String result = zip ? FileHelper.gzFileToString(f) : FileUtils.fileToString(f);
    assertEquals(apples, result);
  }

  /**
   * Test of outStream method, of class NgsResultStreamHandler.
   */
  public void testOutStream() throws IOException {
    NgsResultStreamHandler nr = new NgsResultStreamHandler(mDir, true);
    checkStream(nr.outStream(), NgsResultStreamHandler.OUT_SUFFIX, true);
    nr = new NgsResultStreamHandler(mDir, false);
    checkStream(nr.outStream(), NgsResultStreamHandler.OUT_SUFFIX, false);
  }

  /**
   * Test of samStream method, of class NgsResultStreamHandler.
   */
  public void testSamStream() throws IOException {
    NgsResultStreamHandler nr = new NgsResultStreamHandler(mDir, true);
    checkStream(nr.matedSamStream(), NgsOutputParams.MATED_SAM_FILE_NAME, true);
    nr = new NgsResultStreamHandler(mDir, false);
    checkStream(nr.matedSamStream(), NgsOutputParams.MATED_SAM_FILE_NAME, false);
  }

  /**
   * Test of repeatStream method, of class NgsResultStreamHandler.
   */
  public void testRepeatStream() throws IOException {
    NgsResultStreamHandler nr = new NgsResultStreamHandler(mDir, true);
    checkStream(nr.repeatStream(), NgsOutputParams.REPEATS_FILE_NAME, true);
    nr = new NgsResultStreamHandler(mDir, false);
    checkStream(nr.repeatStream(), NgsOutputParams.REPEATS_FILE_NAME, false);

  }

  /**
   * Test of unmappedStream method, of class NgsResultStreamHandler.
   */
  public void testUnmappedStream() throws IOException {
    NgsResultStreamHandler nr = new NgsResultStreamHandler(mDir, true);
    checkStream(nr.unmappedStream(), NgsOutputParams.UNMAPPED_FILE_NAME, true);
    nr = new NgsResultStreamHandler(mDir, false);
    checkStream(nr.unmappedStream(), NgsOutputParams.UNMAPPED_FILE_NAME, false);

  }

}
