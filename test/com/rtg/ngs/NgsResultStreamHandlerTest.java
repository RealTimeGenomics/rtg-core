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
package com.rtg.ngs;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;

import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;

/**
 * Test class
 */
public class NgsResultStreamHandlerTest extends TestCase {

  public NgsResultStreamHandlerTest(String testName) {
    super(testName);
  }

  File mDir = null;

  @Override
  public void setUp() throws Exception {
    Diagnostic.setLogStream();
    mDir = FileUtils.createTempDir("ngsresultstreamhandler", "test");
    assertTrue(mDir.delete());
  }

  @Override
  public void tearDown() {
    assertTrue(!mDir.exists() || FileHelper.deleteAll(mDir));
    mDir = null;
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
