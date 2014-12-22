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

import com.rtg.util.test.FileHelper;


/**
 */
public class GzipAsynchOutputStreamTest extends AsynchOutputStreamTest {

  public void testNullFile() throws IOException {
    try {
      GzipAsynchOutputStream out = new GzipAsynchOutputStream(null);
      try {
        fail("IllegalArgumentException expected");
      } finally {
        out.close();
      }
    } catch (IllegalArgumentException e) {
      assertEquals("File cannot be null", e.getMessage());
    }
  }

  @Override
  public void testFlush() throws IOException {
    final File file = File.createTempFile("test", "gzipasynch");
    try {
      GzipAsynchOutputStream out = new GzipAsynchOutputStream(file, 1024, 1024);
      try {
        for (int i = 0; i < 1028; i++) {
          out.write((int) 'a');
        }
        out.flush();
        assertEquals(0, out.mQueue.available());
        out.write((int) 'b');
      } finally {
        out.close();
      }
      final String contents = FileHelper.gzFileToString(file);
      assertTrue(contents.startsWith("aaaaaaa"));
      assertTrue(contents.endsWith("aaab"));
      assertEquals(1028 + 1, contents.length());
    } finally {
      assertTrue(FileHelper.deleteAll(file));
    }
  }
}
